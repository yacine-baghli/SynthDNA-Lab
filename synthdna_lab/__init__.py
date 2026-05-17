"""
SynthDNA Lab — Multi-Vendor Physically-Grounded Synthetic DNA Generator.

A production-ready pipeline for generating ultra-realistic synthetic DNA
datasets calibrated to real-world synthesis vendor error profiles.

Supported technology classes:
  - Chemical Phosphoramidite (Twist, IDT, GenScript)
  - Enzymatic TdT (DNA Script, Ansa Biotechnologies)
  - Photolithographic (ETH Zurich Grass/Stark lab)
  - Thermal Chip (Evonetix)

Supported sequencing backends:
  - Illumina SBS
  - Oxford Nanopore (with homopolymer collapse model)
  - PacBio HiFi
"""

__version__ = "1.0.0"
__author__ = "Yacine Baghli"

from .config import PipelineConfig, CompanyProfile, SequencingBackend
from .config import COMPANY_PROFILES, SEQUENCING_BACKENDS, TECH_CLASSES
from .pipeline import TwistRealisticGenerator, generate_dataset
