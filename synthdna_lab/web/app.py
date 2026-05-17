"""
Flask Web Application for SynthDNA Lab.
Premium UI for configuring and launching DNA dataset generation.
"""

import os
import json
import uuid
import threading
import time
from flask import Flask, render_template, request, jsonify, send_file

from synthdna_lab.config import (
    PipelineConfig, COMPANY_PROFILES, SEQUENCING_BACKENDS, TECH_CLASSES
)
from synthdna_lab.pipeline import TwistRealisticGenerator, generate_dataset

_WEB_DIR = os.path.dirname(os.path.abspath(__file__))

app = Flask(__name__,
            template_folder=os.path.join(_WEB_DIR, 'templates'),
            static_folder=os.path.join(_WEB_DIR, 'static'))

# Job tracking
jobs = {}
OUTPUT_DIR = os.path.join(os.getcwd(), 'generated_datasets')
os.makedirs(OUTPUT_DIR, exist_ok=True)


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/api/profiles', methods=['GET'])
def get_profiles():
    profiles = {}
    for key, profile in COMPANY_PROFILES.items():
        profiles[key] = profile.to_dict()
    return jsonify({'profiles': profiles, 'tech_classes': TECH_CLASSES})


@app.route('/api/sequencing_backends', methods=['GET'])
def get_sequencing_backends():
    backends = {}
    for key, backend in SEQUENCING_BACKENDS.items():
        backends[key] = backend.to_dict()
    return jsonify(backends)


@app.route('/api/preview', methods=['POST'])
def preview_sample():
    data = request.json
    profile_key = data.get('profile', 'twist')
    seq_backend = data.get('sequencing_backend', 'illumina')
    target_len = int(data.get('target_len', 110))

    config = PipelineConfig(target_len=target_len)
    config.set_profile(profile_key)
    config.set_sequencing(seq_backend)
    gen = TwistRealisticGenerator(config)

    seed = int(time.time() * 1000) % (2**31)
    center, traces = gen.generate_sample(index=seed, seed_offset=42,
                                          min_traces=3, max_traces=8)

    gc = (center.count('G') + center.count('C')) / len(center)
    max_hp = 1
    cur = 1
    for i in range(1, len(center)):
        if center[i] == center[i-1]:
            cur += 1
            max_hp = max(max_hp, cur)
        else:
            cur = 1

    trace_stats = []
    for t in traces:
        t_gc = (t.count('G') + t.count('C')) / max(len(t), 1)
        len_diff = len(center) - len(t)
        trace_stats.append({
            'sequence': t, 'length': len(t),
            'gc_content': round(t_gc, 4), 'len_diff': len_diff,
        })

    return jsonify({
        'center': center, 'center_len': len(center),
        'gc_content': round(gc, 4), 'max_homopolymer': max_hp,
        'num_traces': len(traces), 'traces': trace_stats,
        'profile': profile_key, 'sequencing_backend': seq_backend,
    })


@app.route('/api/compare', methods=['POST'])
def compare_profiles():
    data = request.json
    profile_a = data.get('profile_a', 'twist')
    profile_b = data.get('profile_b', 'photolitho_ethz')
    target_len = int(data.get('target_len', 110))
    n_samples = min(int(data.get('n_samples', 200)), 500)

    results = {}
    for key in [profile_a, profile_b]:
        config = PipelineConfig(target_len=target_len)
        config.set_profile(key)
        gen = TwistRealisticGenerator(config)

        len_diffs = []
        trace_counts = []
        for i in range(n_samples):
            center, traces = gen.generate_sample(i, seed_offset=99, min_traces=3, max_traces=10)
            trace_counts.append(len(traces))
            for t in traces:
                len_diffs.append(len(center) - len(t))

        p = COMPANY_PROFILES[key]
        total = p.synth_p_del + p.synth_p_sub + p.synth_p_ins
        results[key] = {
            'name': p.name, 'tech_class': p.tech_class,
            'total_error_rate': total,
            'del_frac': p.synth_p_del / total if total > 0 else 0,
            'sub_frac': p.synth_p_sub / total if total > 0 else 0,
            'ins_frac': p.synth_p_ins / total if total > 0 else 0,
            'avg_len_diff': round(sum(len_diffs) / max(len(len_diffs), 1), 2),
            'avg_traces': round(sum(trace_counts) / max(len(trace_counts), 1), 1),
            'gc_range': list(p.gc_range),
            'max_homopolymer': p.max_homopolymer,
            'max_oligo_len': p.max_oligo_len,
        }
    return jsonify(results)


@app.route('/api/generate', methods=['POST'])
def start_generation():
    data = request.json
    profile_key = data.get('profile', 'twist')
    seq_backend = data.get('sequencing_backend', 'illumina')
    target_len = int(data.get('target_len', 110))
    dataset_size = int(data.get('dataset_size', 100000))
    output_format = data.get('output_format', 'fasta')
    min_traces = int(data.get('min_traces', 3))
    max_traces = int(data.get('max_traces', 10))
    num_workers = int(data.get('num_workers', 4))

    job_id = str(uuid.uuid4())[:8]
    jobs[job_id] = {
        'status': 'running', 'progress': 0,
        'message': 'Initializing...', 'result': None, 'error': None,
    }

    def run_job():
        try:
            config = PipelineConfig(
                target_len=target_len, dataset_size=dataset_size,
                min_traces=min_traces, max_traces=max_traces,
                num_workers=num_workers,
            )
            config.set_profile(profile_key)
            config.set_sequencing(seq_backend)

            def progress_cb(pct, msg):
                jobs[job_id]['progress'] = pct
                jobs[job_id]['message'] = msg

            result = generate_dataset(
                config=config, profile_key=profile_key,
                output_dir=OUTPUT_DIR, output_format=output_format,
                progress_callback=progress_cb,
            )
            jobs[job_id]['status'] = 'complete'
            jobs[job_id]['progress'] = 100
            jobs[job_id]['result'] = result
        except Exception as e:
            jobs[job_id]['status'] = 'error'
            jobs[job_id]['error'] = str(e)
            jobs[job_id]['message'] = f"Error: {e}"

    threading.Thread(target=run_job, daemon=True).start()
    return jsonify({'job_id': job_id})


@app.route('/api/status/<job_id>', methods=['GET'])
def job_status(job_id):
    if job_id not in jobs:
        return jsonify({'error': 'Job not found'}), 404
    return jsonify(jobs[job_id])


@app.route('/api/download/<job_id>', methods=['GET'])
def download_result(job_id):
    if job_id not in jobs or jobs[job_id]['status'] != 'complete':
        return jsonify({'error': 'Job not ready'}), 404
    return send_file(jobs[job_id]['result']['output_file'], as_attachment=True)


def main():
    """CLI entry point."""
    print("\n[SynthDNA Lab] Starting web server...")
    print(f"  Profiles: {list(COMPANY_PROFILES.keys())}")
    print(f"  Seq backends: {list(SEQUENCING_BACKENDS.keys())}")
    print(f"  Output dir: {OUTPUT_DIR}")
    print(f"\n  Open http://127.0.0.1:5000\n")
    app.run(debug=True, host='127.0.0.1', port=5000)


if __name__ == '__main__':
    main()
