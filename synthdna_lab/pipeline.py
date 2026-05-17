"""
Pipeline Orchestrator.

Brings together all modules into a single callable generator.
Supports multiprocessing for large dataset generation and exports
to multiple formats (FASTA, CSV, PyTorch .pt).
"""

import random
import math
import time
import os
import json
from typing import Tuple, List, Optional, Dict, Any
from .config import PipelineConfig, COMPANY_PROFILES, SEQUENCING_BACKENDS
from .sequence_generator import TwistSequenceGenerator
from .error_injector import TwistErrorInjector
from .trace_generator import TracePoolGenerator


class TwistRealisticGenerator:
    """
    Production-ready synthetic DNA dataset generator.
    
    Generates (center_sequence, List[corrupted_traces]) pairs that
    faithfully reproduce Twist Bioscience's DNA storage channel.
    """

    def __init__(self, config: Optional[PipelineConfig] = None):
        self.config = config or PipelineConfig()
        self.seq_gen = TwistSequenceGenerator(self.config)
        self.error_inj = TwistErrorInjector(self.config)
        self.pool_gen = TracePoolGenerator(self.config)

    def generate_sample(self, index: int, seed_offset: int = 42,
                        min_traces: int = 3, max_traces: int = 10
                        ) -> Tuple[str, List[str]]:
        """Generate one (center, traces) sample."""
        seed = index + seed_offset

        # Set difficulty from index (curriculum-compatible)
        mode_prob = (index % 100) / 100.0
        if mode_prob < 0.4:
            self.error_inj.set_difficulty(1.0)
        elif mode_prob < 0.7:
            self.error_inj.set_difficulty(3.0)
        elif mode_prob < 0.9:
            self.error_inj.set_difficulty(6.0)
        else:
            self.error_inj.set_difficulty(1.0)

        # Override pool gen limits
        self.pool_gen.min_traces = min_traces
        self.pool_gen.max_traces = max_traces

        center = self.seq_gen.generate(seed)
        traces = self.pool_gen.generate_pool(center, self.error_inj, seed)
        return center, traces


def _generate_chunk(args) -> Dict[str, Any]:
    """Worker function for parallel generation."""
    start_idx, end_idx, config_dict, seed_offset, min_traces, max_traces = args
    
    config = PipelineConfig()
    config.target_len = config_dict['target_len']
    config.set_profile(config_dict['profile_key'])
    config.set_sequencing(config_dict.get('seq_backend', 'illumina'))
    
    gen = TwistRealisticGenerator(config)
    
    centers = []
    clusters = []
    stats = {'total_bases': 0, 'total_deletions': 0, 'total_insertions': 0,
             'total_substitutions': 0, 'total_traces': 0}

    for i in range(start_idx, end_idx):
        center, traces = gen.generate_sample(i, seed_offset, min_traces, max_traces)
        centers.append(center)
        clusters.append(traces)
        stats['total_bases'] += len(center) * len(traces)
        stats['total_traces'] += len(traces)
        # Approximate error counting (length diff = net indels)
        for t in traces:
            diff = len(center) - len(t)
            if diff > 0:
                stats['total_deletions'] += diff
            elif diff < 0:
                stats['total_insertions'] += abs(diff)

    return {'centers': centers, 'clusters': clusters, 'stats': stats}


def generate_dataset(config: PipelineConfig, profile_key: str = "twist",
                     output_dir: str = ".", output_format: str = "fasta",
                     progress_callback=None) -> Dict[str, Any]:
    """
    Generate a full dataset with progress reporting.
    
    Args:
        config: Pipeline configuration
        profile_key: Company profile key
        output_dir: Directory to save output files
        output_format: "fasta", "csv", or "pt"
        progress_callback: Optional callable(percent, message)
    
    Returns:
        Dictionary with generation statistics
    """
    import multiprocessing as mp

    config.set_profile(profile_key)
    n = config.dataset_size
    num_workers = min(config.num_workers, max(1, mp.cpu_count() - 1))
    
    os.makedirs(output_dir, exist_ok=True)
    
    if progress_callback:
        progress_callback(0, f"Starting generation of {n:,} samples with {num_workers} workers...")
    
    start_time = time.time()
    
    chunk_size = math.ceil(n / num_workers)
    tasks = []
    config_dict = {
        'target_len': config.target_len,
        'profile_key': profile_key,
        'seq_backend': config.sequencing.key,
    }
    
    for w in range(num_workers):
        s = w * chunk_size
        e = min((w + 1) * chunk_size, n)
        if s < e:
            tasks.append((s, e, config_dict, config.seed,
                          config.min_traces, config.max_traces))
    
    # Use multiprocessing for large datasets, single-process for small
    all_centers = []
    all_clusters = []
    agg_stats = {'total_bases': 0, 'total_deletions': 0,
                 'total_insertions': 0, 'total_traces': 0}
    
    if n >= 10000 and num_workers > 1:
        with mp.Pool(num_workers) as pool:
            for i, result in enumerate(pool.imap_unordered(_generate_chunk, tasks)):
                all_centers.extend(result['centers'])
                all_clusters.extend(result['clusters'])
                for k in agg_stats:
                    agg_stats[k] += result['stats'][k]
                pct = int((i + 1) / len(tasks) * 80)
                if progress_callback:
                    progress_callback(pct, f"Generated {len(all_centers):,}/{n:,} samples...")
    else:
        for i, task in enumerate(tasks):
            result = _generate_chunk(task)
            all_centers.extend(result['centers'])
            all_clusters.extend(result['clusters'])
            for k in agg_stats:
                agg_stats[k] += result['stats'][k]
            if progress_callback:
                pct = int((i + 1) / len(tasks) * 80)
                progress_callback(pct, f"Generated {len(all_centers):,}/{n:,} samples...")
    
    gen_time = time.time() - start_time
    
    if progress_callback:
        progress_callback(85, "Saving dataset...")
    
    # Save to requested format
    filename = f"synthetic_{profile_key}_{n}"
    
    if output_format == "fasta":
        filepath = os.path.join(output_dir, filename + ".fasta")
        _save_fasta(filepath, all_centers, all_clusters)
    elif output_format == "csv":
        filepath = os.path.join(output_dir, filename + ".csv")
        _save_csv(filepath, all_centers, all_clusters)
    elif output_format == "pt":
        import torch
        filepath = os.path.join(output_dir, filename + ".pt")
        torch.save({'centers': all_centers, 'clusters': all_clusters}, filepath)
    else:
        filepath = os.path.join(output_dir, filename + ".fasta")
        _save_fasta(filepath, all_centers, all_clusters)
    
    file_size = os.path.getsize(filepath)
    
    if progress_callback:
        progress_callback(100, "Complete!")
    
    # Compute summary statistics
    gc_values = []
    trace_counts = []
    for c in all_centers[:1000]:  # Sample first 1000 for speed
        gc = (c.count('G') + c.count('C')) / len(c)
        gc_values.append(gc)
    for cl in all_clusters[:1000]:
        trace_counts.append(len(cl))
    
    summary = {
        'profile': profile_key,
        'profile_name': config.profile.name,
        'dataset_size': len(all_centers),
        'target_len': config.target_len,
        'generation_time_sec': round(gen_time, 2),
        'samples_per_sec': round(len(all_centers) / gen_time, 1),
        'output_file': filepath,
        'output_format': output_format,
        'file_size_mb': round(file_size / (1024 * 1024), 2),
        'avg_gc_content': round(sum(gc_values) / max(len(gc_values), 1), 4),
        'avg_traces_per_sample': round(sum(trace_counts) / max(len(trace_counts), 1), 2),
        'total_traces_generated': agg_stats['total_traces'],
        'approx_deletion_rate': round(agg_stats['total_deletions'] / max(agg_stats['total_bases'], 1), 6),
        'approx_insertion_rate': round(agg_stats['total_insertions'] / max(agg_stats['total_bases'], 1), 6),
    }
    
    # Save metadata
    meta_path = os.path.join(output_dir, filename + "_metadata.json")
    with open(meta_path, 'w') as f:
        json.dump(summary, f, indent=2)
    
    return summary


def _save_fasta(path: str, centers: List[str], clusters: List[List[str]]):
    """Save dataset in FASTA format."""
    with open(path, 'w') as f:
        for i, (center, traces) in enumerate(zip(centers, clusters)):
            f.write(f">sample_{i}_center\n{center}\n")
            for j, trace in enumerate(traces):
                f.write(f">sample_{i}_trace_{j}\n{trace}\n")


def _save_csv(path: str, centers: List[str], clusters: List[List[str]]):
    """Save dataset in CSV format."""
    with open(path, 'w') as f:
        f.write("sample_id,type,index,sequence\n")
        for i, (center, traces) in enumerate(zip(centers, clusters)):
            f.write(f"{i},center,0,{center}\n")
            for j, trace in enumerate(traces):
                f.write(f"{i},trace,{j},{trace}\n")
