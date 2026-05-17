"""
Synthesis Company Profiles, Sequencing Backends & Pipeline Configuration.

Each CompanyProfile encapsulates empirically-calibrated error rates,
sequence constraints, and coverage statistics for a specific DNA
synthesis vendor's platform.

Each SequencingBackend captures the error signature of a specific
sequencing technology (Illumina SBS, Oxford Nanopore, PacBio HiFi).

Profiles are organized by technology class:
  - chemical:      Phosphoramidite column/array (Twist, IDT, etc.)
  - enzymatic:     TdT-based synthesis (DNA Script, Ansa, Telesis Bio)
  - photolitho:    Photolithographic array (ETH Zurich Grass/Stark)
  - thermal_chip:  Thermal MEMS (Evonetix)
"""

from dataclasses import dataclass, field, asdict
from typing import Tuple, Dict, Any, List
import json


# ═══════════════════════════════════════════════════════════════
# Sequencing Backend
# ═══════════════════════════════════════════════════════════════

@dataclass
class SequencingBackend:
    """A sequencing technology's error signature."""
    name: str = "Illumina SBS"
    key: str = "illumina"
    description: str = "Illumina sequencing-by-synthesis. Substitution-dominated, negligible indels."
    seq_p_sub: float = 1.0e-3
    seq_p_del: float = 0.0
    seq_p_ins: float = 0.0
    seq_position_decay: float = 0.5
    seq_transition_bias: float = 0.6
    # Nanopore-specific: systematic homopolymer length miscall
    homopolymer_collapse: bool = False
    homopolymer_miscall_rate: float = 0.0  # P(miscall) for runs >= 4bp

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


SEQUENCING_BACKENDS: Dict[str, SequencingBackend] = {
    "illumina": SequencingBackend(
        name="Illumina SBS",
        key="illumina",
        description=(
            "Illumina sequencing-by-synthesis with fluorescent reversible terminators. "
            "Substitution-dominated (~0.1% per base), negligible indels. "
            "Quality decays toward 3' end. Transitions (A<>G, C<>T) enriched."
        ),
        seq_p_sub=1.0e-3,
        seq_p_del=0.0,
        seq_p_ins=0.0,
        seq_position_decay=0.5,
        seq_transition_bias=0.6,
    ),
    "nanopore": SequencingBackend(
        name="Oxford Nanopore",
        key="nanopore",
        description=(
            "Oxford Nanopore MinION/PromethION. Indel-dominated errors (~3-5% total). "
            "Systematic homopolymer length miscalls: runs >= 4bp often under-reported "
            "by 1-2 bases. Position-independent. Critical for DNA storage read simulation."
        ),
        seq_p_sub=5.0e-3,
        seq_p_del=2.0e-2,
        seq_p_ins=1.5e-2,
        seq_position_decay=0.1,
        seq_transition_bias=0.5,
        homopolymer_collapse=True,
        homopolymer_miscall_rate=0.25,  # 25% chance of +-1 miscall for runs >= 4bp
    ),
    "pacbio": SequencingBackend(
        name="PacBio HiFi",
        key="pacbio",
        description=(
            "PacBio HiFi (CCS) sequencing. Ultra-low error rate (~0.01-0.03%). "
            "Balanced error types. Minimal position dependence. "
            "Gold standard for long-read accuracy."
        ),
        seq_p_sub=1.0e-4,
        seq_p_del=5.0e-5,
        seq_p_ins=5.0e-5,
        seq_position_decay=0.05,
        seq_transition_bias=0.5,
    ),
}


# ═══════════════════════════════════════════════════════════════
# Company Profile
# ═══════════════════════════════════════════════════════════════

@dataclass
class CompanyProfile:
    """A vendor-specific DNA synthesis error & constraint profile."""
    name: str = "Twist Bioscience"
    platform: str = "Silicon-based phosphoramidite"
    tech_class: str = "chemical"  # chemical | enzymatic | photolitho | thermal_chip
    description: str = ""

    # ── Sequence Constraints ──
    gc_range: Tuple[float, float] = (0.25, 0.75)
    gc_optimal: Tuple[float, float] = (0.35, 0.65)
    max_homopolymer: int = 6
    max_dinuc_repeat: int = 10
    min_oligo_len: int = 20
    max_oligo_len: int = 300

    # ── Synthesis Errors (Layer 1 — dominant) ──
    synth_p_del: float = 2.5e-4
    synth_p_sub: float = 1.2e-4
    synth_p_ins: float = 0.4e-4
    synth_position_slope: float = 0.3
    homopoly_alpha: float = 0.6
    homopoly_beta: float = 0.35

    # ── PCR Amplification Errors (Layer 2) ──
    pcr_error_per_base: float = 5.0e-6
    pcr_sub_frac: float = 0.85
    pcr_del_frac: float = 0.10
    pcr_ins_frac: float = 0.05
    pcr_cycles: int = 18

    # ── Sequencing Errors (Layer 3 — defaults, overridden by backend) ──
    seq_p_sub: float = 1.0e-3
    seq_position_decay: float = 0.5
    seq_transition_bias: float = 0.6

    # ── Coverage Distribution ──
    mean_coverage: float = 8.0
    nb_dispersion: float = 2.0
    gc_bias_strength: float = 0.05

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        # Tuples become lists in asdict; keep them as lists for JSON
        return d

    def to_json(self) -> str:
        return json.dumps(self.to_dict(), indent=2)


# ═══════════════════════════════════════════════════════════════
# Pre-built vendor profiles
# ═══════════════════════════════════════════════════════════════

COMPANY_PROFILES: Dict[str, CompanyProfile] = {

    # ─── Chemical Phosphoramidite: Silicon Array ───
    "twist": CompanyProfile(
        name="Twist Bioscience",
        platform="Silicon-based phosphoramidite array",
        tech_class="chemical",
        description=(
            "Twist Bioscience silicon-based high-throughput oligo synthesis. "
            "Industry-leading error rate ~1:2000-1:3000 nt. Deletions dominate "
            "(~60% of errors) due to coupling failure. PCR with KAPA HiFi, "
            "sequenced on Illumina SBS. The gold standard for DNA data storage."
        ),
    ),

    # ─── Chemical Phosphoramidite: Column ───
    "idt": CompanyProfile(
        name="IDT (Integrated DNA Technologies)",
        platform="Column-based phosphoramidite",
        tech_class="chemical",
        description=(
            "IDT column-based synthesis. Higher per-oligo fidelity than array platforms "
            "(~1:5000 nt for Ultramers). Lower throughput but gold-standard purity. "
            "Owned by Danaher. Used extensively in CRISPR and diagnostics."
        ),
        synth_p_del=1.5e-4,
        synth_p_sub=0.8e-4,
        synth_p_ins=0.3e-4,
        synth_position_slope=0.2,
        max_oligo_len=200,
        max_homopolymer=8,
    ),

    # ─── Chemical Phosphoramidite: CMOS Array ───
    "genscript": CompanyProfile(
        name="GenScript Biotech",
        platform="CMOS array-based phosphoramidite",
        tech_class="chemical",
        description=(
            "GenScript semiconductor chip-based oligo pool synthesis. "
            "Array platform with error rate ~1:200-1:500 per base for pools. "
            "100% sequence-verified for clonal gene synthesis. "
            "Massive parallelism for CRISPR libraries."
        ),
        synth_p_del=3.0e-3,
        synth_p_sub=1.5e-3,
        synth_p_ins=0.5e-3,
        synth_position_slope=0.35,
        homopoly_alpha=0.7,
        max_oligo_len=200,
        max_homopolymer=6,
        gc_range=(0.20, 0.80),
    ),

    # ─── Enzymatic TdT: On-demand benchtop ───
    "dnascript": CompanyProfile(
        name="DNA Script (SYNTAX)",
        platform="Enzymatic TdT synthesis",
        tech_class="enzymatic",
        description=(
            "DNA Script SYNTAX platform. Template-independent enzymatic synthesis "
            "using engineered TdT with reversibly terminated nucleotides. "
            "Aqueous chemistry - handles poly(A), high GC, ITRs. ~23% perfect clones. "
            "Error profile shifted toward balanced del/sub with 2x higher insertions."
        ),
        gc_range=(0.10, 0.90),
        gc_optimal=(0.30, 0.70),
        max_homopolymer=15,
        max_dinuc_repeat=20,
        min_oligo_len=20,
        max_oligo_len=500,
        synth_p_del=4.0e-3,
        synth_p_sub=4.0e-3,
        synth_p_ins=2.0e-3,
        synth_position_slope=0.15,
        homopoly_alpha=0.3,
        homopoly_beta=0.20,
        gc_bias_strength=0.03,
    ),

    # ─── Enzymatic TdT: High-fidelity clonal ───
    "ansa": CompanyProfile(
        name="Ansa Biotechnologies",
        platform="Enzymatic tethered-TdT",
        tech_class="enzymatic",
        description=(
            "Ansa Bio enzymatic synthesis with tethered-linker TdT strategy. "
            "Stepwise yield >99.9%. Error rate ~1:750-1:1000. "
            "Handles extreme GC (10-94%), hairpins, repeats. "
            "Clonal DNA up to 50 kb. Long-read QC verified."
        ),
        gc_range=(0.10, 0.94),
        gc_optimal=(0.30, 0.70),
        max_homopolymer=12,
        max_dinuc_repeat=15,
        max_oligo_len=1000,
        synth_p_del=5.0e-4,
        synth_p_sub=5.0e-4,
        synth_p_ins=3.0e-4,
        synth_position_slope=0.10,
        homopoly_alpha=0.25,
        homopoly_beta=0.18,
        gc_bias_strength=0.02,
    ),

    # ─── Photolithographic Array ───
    "photolitho_ethz": CompanyProfile(
        name="ETH Zurich Photolithographic",
        platform="Photolithographic UV array",
        tech_class="photolitho",
        description=(
            "Photolithographic array synthesis (Grass/Stark/Heckel labs, ETH Zurich 2024). "
            "Ultra-low-cost but ~12% error per nt. Deletion-dominated. "
            "Only ~2% of raw reads are error-free. Requires sophisticated ECC. "
            "THE critical frontier for DNA data storage cost reduction."
        ),
        gc_range=(0.15, 0.85),
        gc_optimal=(0.30, 0.70),
        max_homopolymer=4,
        max_dinuc_repeat=6,
        max_oligo_len=150,
        synth_p_del=0.082,
        synth_p_sub=0.025,
        synth_p_ins=0.016,
        synth_position_slope=0.5,
        homopoly_alpha=1.2,
        homopoly_beta=0.45,
        nb_dispersion=1.5,
        gc_bias_strength=0.10,
    ),

    # ─── Thermal Chip ───
    "evonetix": CompanyProfile(
        name="Evonetix",
        platform="Thermal MEMS chip + Binary Assembly",
        tech_class="thermal_chip",
        description=(
            "Evonetix semiconductor chip with independent thermal control per pixel. "
            "Binary Assembly removes errors via mismatch melting temperature discrimination. "
            "ML-predicted Tm for error rejection. Claims 'orders of magnitude better' than "
            "conventional. Estimated <1e-5 per nt error rate."
        ),
        gc_range=(0.20, 0.80),
        gc_optimal=(0.35, 0.65),
        max_homopolymer=10,
        max_dinuc_repeat=12,
        max_oligo_len=1000,
        synth_p_del=5.0e-6,
        synth_p_sub=3.0e-6,
        synth_p_ins=2.0e-6,
        synth_position_slope=0.05,
        homopoly_alpha=0.2,
        homopoly_beta=0.15,
        gc_bias_strength=0.02,
    ),
}

# Technology class metadata for UI rendering
TECH_CLASSES = {
    "chemical": {"label": "Chemical Phosphoramidite", "color": "#06d6a0", "icon": "flask"},
    "enzymatic": {"label": "Enzymatic TdT", "color": "#8338ec", "icon": "dna"},
    "photolitho": {"label": "Photolithographic", "color": "#ef476f", "icon": "sun"},
    "thermal_chip": {"label": "Thermal Chip", "color": "#ffd166", "icon": "cpu"},
}


# ═══════════════════════════════════════════════════════════════
# Pipeline Configuration
# ═══════════════════════════════════════════════════════════════

@dataclass
class PipelineConfig:
    """Full pipeline configuration combining profile + user settings."""
    target_len: int = 110
    dataset_size: int = 100_000
    min_traces: int = 3
    max_traces: int = 10
    markov_order: int = 3
    num_workers: int = 4
    seed: int = 42

    profile: CompanyProfile = field(default_factory=lambda: COMPANY_PROFILES["twist"])
    sequencing: SequencingBackend = field(default_factory=lambda: SEQUENCING_BACKENDS["illumina"])

    def set_profile(self, profile_key: str):
        if profile_key in COMPANY_PROFILES:
            self.profile = COMPANY_PROFILES[profile_key]
        else:
            raise ValueError(f"Unknown profile: {profile_key}. Available: {list(COMPANY_PROFILES.keys())}")

    def set_sequencing(self, backend_key: str):
        if backend_key in SEQUENCING_BACKENDS:
            self.sequencing = SEQUENCING_BACKENDS[backend_key]
        else:
            raise ValueError(f"Unknown backend: {backend_key}. Available: {list(SEQUENCING_BACKENDS.keys())}")
