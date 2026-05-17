"""
Vercel Serverless Entry Point for SynthDNA Lab.

This module wraps the Flask application as a Vercel serverless function.
Vercel calls the `app` WSGI application for all incoming requests.
"""

import sys
import os

# Add project root to Python path so synthdna_lab is importable
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import json
import time
from flask import Flask, render_template, request, jsonify

from synthdna_lab.config import (
    PipelineConfig, COMPANY_PROFILES, SEQUENCING_BACKENDS, TECH_CLASSES
)
from synthdna_lab.pipeline import TwistRealisticGenerator

# ── Flask app ──
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_TEMPLATE_DIR = os.path.join(_ROOT, 'synthdna_lab', 'web', 'templates')
_STATIC_DIR = os.path.join(_ROOT, 'public', 'static')

app = Flask(__name__,
            template_folder=_TEMPLATE_DIR,
            static_folder=_STATIC_DIR,
            static_url_path='/static')


# ═══════════════════════════════════════
# Routes
# ═══════════════════════════════════════

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
    data = request.json or {}
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
    data = request.json or {}
    profile_a = data.get('profile_a', 'twist')
    profile_b = data.get('profile_b', 'photolitho_ethz')
    target_len = int(data.get('target_len', 110))
    # Reduced sample count for serverless (timeout constraint)
    n_samples = min(int(data.get('n_samples', 100)), 200)

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
    """
    Serverless generation — returns samples directly in JSON.
    No background jobs or filesystem writes on Vercel.
    """
    data = request.json or {}
    profile_key = data.get('profile', 'twist')
    seq_backend = data.get('sequencing_backend', 'illumina')
    target_len = int(data.get('target_len', 110))
    # Hard cap for serverless
    dataset_size = min(int(data.get('dataset_size', 100)), 500)

    config = PipelineConfig(target_len=target_len)
    config.set_profile(profile_key)
    config.set_sequencing(seq_backend)
    gen = TwistRealisticGenerator(config)

    start = time.time()
    centers = []
    clusters = []
    for i in range(dataset_size):
        center, traces = gen.generate_sample(i, seed_offset=42, min_traces=3, max_traces=10)
        centers.append(center)
        clusters.append(traces)

    elapsed = time.time() - start

    # Build FASTA content as downloadable string
    fasta_lines = []
    for i, (c, ts) in enumerate(zip(centers, clusters)):
        fasta_lines.append(f">sample_{i}_center\n{c}")
        for j, t in enumerate(ts):
            fasta_lines.append(f">sample_{i}_trace_{j}\n{t}")

    return jsonify({
        'status': 'complete',
        'dataset_size': dataset_size,
        'generation_time_sec': round(elapsed, 2),
        'samples_per_sec': round(dataset_size / max(elapsed, 0.01), 1),
        'fasta_content': '\n'.join(fasta_lines),
        'note': 'Serverless mode: max 500 samples. Use local CLI for larger datasets.',
    })
