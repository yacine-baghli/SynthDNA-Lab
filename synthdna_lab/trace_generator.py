"""
Module 2: Multi-Copy Trace Generator.

Simulates a PCR amplification pool with Negative Binomial coverage
distribution per source sequence. Each trace is independently
corrupted by the TwistErrorInjector.
"""

import random
import numpy as np
from typing import List
from .config import PipelineConfig
from .error_injector import TwistErrorInjector


class TracePoolGenerator:
    """
    Simulates a PCR amplification pool with realistic coverage
    distribution per source sequence.
    
    Coverage is modeled as NegBin(r, p) where:
      - mean = mean_coverage * gc_penalty
      - variance = mean + mean²/r  (overdispersed Poisson)
    """

    def __init__(self, config: PipelineConfig):
        self.min_traces = config.min_traces
        self.max_traces = config.max_traces
        self.mean_coverage = config.profile.mean_coverage
        self.dispersion = config.profile.nb_dispersion
        self.gc_bias_strength = config.profile.gc_bias_strength

    def sample_coverage(self, gc_content: float, seed: int) -> int:
        """
        Sample number of traces from a GC-modulated Negative Binomial.
        
        The NB is parameterized as:
          r = dispersion (shape parameter)
          p = r / (r + μ_adj)
          Then n ~ NegBin(r, p) has mean = μ_adj
        """
        rng = np.random.RandomState(seed % (2**31))
        
        # GC deviation penalty (mild: sequences near 50% → 1.0x, extreme → 0.9x)
        gc_penalty = 1.0 - self.gc_bias_strength * abs(gc_content - 0.5) * 4.0
        gc_penalty = max(gc_penalty, 0.5)

        mu_adj = self.mean_coverage * gc_penalty
        r = self.dispersion

        # NB parameterization: p = r/(r+mu), sample from NB(r, p)
        p = r / (r + mu_adj)
        n = rng.negative_binomial(r, p)

        return int(np.clip(n, self.min_traces, self.max_traces))

    def generate_pool(self, center: str, error_injector: TwistErrorInjector,
                      seed: int) -> List[str]:
        """
        Generate a pool of corrupted traces for one source sequence.
        Each trace is independently corrupted.
        """
        gc = (center.count('G') + center.count('C')) / max(len(center), 1)
        n_traces = self.sample_coverage(gc, seed)

        traces = []
        for i in range(n_traces):
            # Each trace gets a unique, deterministic RNG
            trace_rng = random.Random(seed + i * 7919 + 1)
            trace = error_injector.corrupt(center, trace_rng)
            traces.append(trace)

        return traces
