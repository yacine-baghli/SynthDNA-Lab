"""
Module 3: Multi-Vendor Error Injector.

Three-layer error injection with sequencing backend selection:
  Layer 1: Vendor-specific synthesis errors (profile-driven)
  Layer 2: PCR amplification (rare HiFi polymerase errors)
  Layer 3: Sequencing backend (Illumina SBS / Nanopore / PacBio HiFi)

Nanopore Layer 3 includes a specialized homopolymer collapse model:
systematic length miscalls where runs >= 4bp are under-reported by 1-2 bases.
This is fundamentally different from random per-base deletions.
"""

import random
import math
import numpy as np
from typing import List
from .config import PipelineConfig, CompanyProfile, SequencingBackend

BASES = ['A', 'C', 'G', 'T']
BASE_TO_IDX = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


class TwistErrorInjector:
    """
    Core error engine with three independent corruption layers.
    Supports all vendor profiles and sequencing backends.
    """

    def __init__(self, config: PipelineConfig):
        self.profile = config.profile
        self.sequencing = config.sequencing
        self.difficulty_scale = 1.0

        self.synth_sub_matrix = self._build_synth_sub_matrix()
        self.seq_sub_matrix = self._build_seq_sub_matrix()

    def _build_synth_sub_matrix(self) -> np.ndarray:
        """Synthesis substitution matrix (depurination bias for chemical, uniform for enzymatic)."""
        mat = np.ones((4, 4), dtype=np.float64)
        if self.profile.tech_class == "chemical":
            mat[0][2] = 1.5  # A->G (depurination)
            mat[2][0] = 1.5  # G->A
        elif self.profile.tech_class == "enzymatic":
            pass  # Uniform - enzymatic errors are less biased
        elif self.profile.tech_class == "photolitho":
            mat[0][2] = 1.3  # Slight purine bias
            mat[2][0] = 1.3
        np.fill_diagonal(mat, 0.0)
        for i in range(4):
            s = mat[i].sum()
            if s > 0:
                mat[i] /= s
        return mat

    def _build_seq_sub_matrix(self) -> np.ndarray:
        """Sequencing substitution matrix (technology-dependent)."""
        bias = self.sequencing.seq_transition_bias
        mat = np.ones((4, 4), dtype=np.float64)
        if bias > 0.5:
            r = bias / (1 - bias)
            mat[0][2] = r  # A->G
            mat[2][0] = r  # G->A
            mat[1][3] = r  # C->T
            mat[3][1] = r  # T->C
        np.fill_diagonal(mat, 0.0)
        for i in range(4):
            s = mat[i].sum()
            if s > 0:
                mat[i] /= s
        return mat

    def _compute_homopolymer_runs(self, seq: str) -> List[int]:
        """For each position, compute the length of its homopolymer run."""
        L = len(seq)
        if L == 0:
            return []
        runs = [1] * L
        run_len = 1
        start = 0
        for j in range(1, L):
            if seq[j] == seq[j - 1]:
                run_len += 1
            else:
                for k in range(start, j):
                    runs[k] = run_len
                run_len = 1
                start = j
        for k in range(start, L):
            runs[k] = run_len
        return runs

    def _position_scale(self, pos: int, seq_len: int) -> float:
        """Position-dependent error scaling for synthesis."""
        slope = self.profile.synth_position_slope
        return 1.0 + slope * (1.0 - pos / max(seq_len - 1, 1))

    def _apply_synthesis_errors(self, seq: str, rng: random.Random) -> str:
        """Layer 1: Vendor-specific synthesis corruption."""
        p = self.profile
        alpha = p.homopoly_alpha
        beta = p.homopoly_beta
        hp_runs = self._compute_homopolymer_runs(seq)
        L = len(seq)
        result = []
        i = 0

        while i < L:
            base = seq[i]
            pos_scale = self._position_scale(i, L)
            hp_len = hp_runs[i]
            hp_scale = 1.0 + alpha * (math.exp(beta * max(hp_len - 2, 0)) - 1.0)

            p_del = min(p.synth_p_del * pos_scale * hp_scale * self.difficulty_scale, 0.5)
            p_ins = min(p.synth_p_ins * pos_scale * self.difficulty_scale, 0.3)
            p_sub = min(p.synth_p_sub * pos_scale * self.difficulty_scale, 0.3)

            r = rng.random()
            if r < p_del:
                i += 1
            elif r < p_del + p_ins:
                result.append(BASES[rng.randint(0, 3)])
            elif r < p_del + p_ins + p_sub:
                bi = BASE_TO_IDX[base]
                result.append(BASES[self._weighted_choice(self.synth_sub_matrix[bi], rng)])
                i += 1
            else:
                result.append(base)
                i += 1

        return ''.join(result)

    def _apply_pcr_errors(self, seq: str, rng: random.Random) -> str:
        """Layer 2: PCR polymerase errors (very rare, mostly substitutions)."""
        p = self.profile
        p_err = p.pcr_error_per_base * self.difficulty_scale
        result = []
        i = 0

        while i < len(seq):
            base = seq[i]
            if rng.random() < p_err:
                r2 = rng.random()
                if r2 < p.pcr_sub_frac:
                    choices = [b for b in BASES if b != base]
                    result.append(rng.choice(choices))
                    i += 1
                elif r2 < p.pcr_sub_frac + p.pcr_del_frac:
                    i += 1
                else:
                    result.append(BASES[rng.randint(0, 3)])
            else:
                result.append(base)
                i += 1

        return ''.join(result)

    def _apply_sequencing_errors(self, seq: str, rng: random.Random) -> str:
        """
        Layer 3: Sequencing backend errors.
        
        Illumina: substitution-only with 3' quality decay.
        Nanopore:  IDS errors + systematic homopolymer length miscalls.
        PacBio:    Very low balanced IDS errors.
        """
        sb = self.sequencing

        # For Nanopore: first apply homopolymer collapse, then per-base errors
        if sb.homopolymer_collapse:
            seq = self._apply_nanopore_homopolymer_collapse(seq, rng)

        L = len(seq)
        result = []
        i = 0

        while i < L:
            base = seq[i]
            # Position-dependent quality decay
            q_scale = 1.0 + sb.seq_position_decay * (i / max(L - 1, 1))
            scale = q_scale * self.difficulty_scale

            p_del = min(sb.seq_p_del * scale, 0.3)
            p_ins = min(sb.seq_p_ins * scale, 0.3)
            p_sub = min(sb.seq_p_sub * scale, 0.3)

            r = rng.random()
            if r < p_del:
                # Deletion
                i += 1
            elif r < p_del + p_ins:
                # Insertion
                result.append(BASES[rng.randint(0, 3)])
            elif r < p_del + p_ins + p_sub:
                # Substitution with bias matrix
                bi = BASE_TO_IDX[base]
                result.append(BASES[self._weighted_choice(self.seq_sub_matrix[bi], rng)])
                i += 1
            else:
                result.append(base)
                i += 1

        return ''.join(result)

    def _apply_nanopore_homopolymer_collapse(self, seq: str, rng: random.Random) -> str:
        """
        Nanopore-specific: Systematic homopolymer length miscalls.
        
        Unlike random per-base deletions, Nanopore's signal-level basecalling
        systematically miscalls homopolymer run lengths. For runs >= 4bp:
          - P(miscall) = homopolymer_miscall_rate (default 25%)
          - When miscall occurs: 70% under-report by 1, 20% under by 2, 10% over by 1
        
        This creates a fundamentally different error signature from synthesis deletions.
        """
        if not seq:
            return seq

        sb = self.sequencing
        result = []
        i = 0
        L = len(seq)

        while i < L:
            # Find homopolymer run starting at i
            base = seq[i]
            run_end = i + 1
            while run_end < L and seq[run_end] == base:
                run_end += 1
            run_len = run_end - i

            if run_len >= 4 and rng.random() < sb.homopolymer_miscall_rate:
                # Systematic miscall
                r = rng.random()
                if r < 0.70:
                    # Under-report by 1 (most common)
                    reported_len = max(run_len - 1, 1)
                elif r < 0.90:
                    # Under-report by 2
                    reported_len = max(run_len - 2, 1)
                else:
                    # Over-report by 1
                    reported_len = run_len + 1
                result.extend([base] * reported_len)
            else:
                # Faithful copy of the run
                result.extend([base] * run_len)

            i = run_end

        return ''.join(result)

    def corrupt(self, center: str, rng: random.Random) -> str:
        """Full pipeline corruption: synthesis -> PCR -> sequencing."""
        seq = self._apply_synthesis_errors(center, rng)
        seq = self._apply_pcr_errors(seq, rng)
        seq = self._apply_sequencing_errors(seq, rng)
        return seq

    def set_difficulty(self, scale: float):
        """Set curriculum difficulty multiplier."""
        self.difficulty_scale = max(0.1, scale)

    def update_config(self, config: PipelineConfig):
        """Re-initialize from updated config (profile or sequencing change)."""
        self.profile = config.profile
        self.sequencing = config.sequencing
        self.synth_sub_matrix = self._build_synth_sub_matrix()
        self.seq_sub_matrix = self._build_seq_sub_matrix()

    @staticmethod
    def _weighted_choice(probs: np.ndarray, rng: random.Random) -> int:
        r = rng.random()
        cumsum = 0.0
        for i in range(len(probs)):
            cumsum += probs[i]
            if r < cumsum:
                return i
        return len(probs) - 1
