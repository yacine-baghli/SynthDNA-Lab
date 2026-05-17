<div align="center">

# SynthDNA Lab

### Multi-Vendor Physically-Grounded Synthetic DNA Dataset Generator

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Flask](https://img.shields.io/badge/Flask-3.0+-black.svg)](https://flask.palletsprojects.com/)

*Generate ultra-realistic synthetic DNA datasets calibrated to real-world synthesis vendor error profiles for DNA data storage research.*

</div>

---

## Overview

**SynthDNA Lab** is a production-ready pipeline for generating synthetic DNA sequence datasets that faithfully reproduce the physical, chemical, and biological error profiles of modern DNA synthesis platforms. Unlike simple random ACGT generation or uniform noise models, SynthDNA Lab implements a **3-layer error injection architecture** calibrated from empirical data published by leading research groups.


<img width="1869" height="1064" alt="image" src="https://github.com/user-attachments/assets/3e2c0515-c4fc-43ee-857a-ab8728d40e61" />


### Why does this matter?

Training machine learning models for **DNA trace reconstruction** (e.g., diffusion models, transformers) requires datasets that capture the true noise characteristics of the DNA data storage channel. SynthDNA Lab provides this by modeling:

1. **Synthesis errors** — vendor-specific, position-dependent, homopolymer-biased
2. **PCR amplification errors** — rare polymerase mistakes during copying
3. **Sequencing errors** — technology-specific (Illumina, Nanopore, PacBio)

---

## Supported Platforms

### Synthesis Vendors (7 profiles across 4 technology classes)

| Technology Class | Vendor | Error Rate | Key Characteristics |
|---|---|---|---|
| **Chemical Phosphoramidite** | Twist Bioscience | ~1:2,500/nt | Silicon array, deletion-dominated, DNA storage gold standard |
| | IDT (Danaher) | ~1:5,000/nt | Column-based, highest per-oligo fidelity |
| | GenScript | ~1:200/nt | CMOS array, oligo pool synthesis |
| **Enzymatic TdT** | DNA Script (SYNTAX) | ~1:100/nt | Aqueous chemistry, handles poly(A)/ITRs |
| | Ansa Biotechnologies | ~1:800/nt | Tethered-TdT, extreme GC tolerance (10-94%) |
| **Photolithographic** | ETH Zurich (Grass/Stark) | ~1:8/nt | Ultra-low-cost, ~2% error-free reads, critical for storage cost reduction |
| **Thermal Chip** | Evonetix | ~1:100,000/nt | Binary Assembly error removal, next-gen accuracy |

### Sequencing Backends (3 technologies)

| Backend | Error Type | Special Features |
|---|---|---|
| **Illumina SBS** | Substitution-dominated (0.1%) | 3' quality decay, transition bias |
| **Oxford Nanopore** | Indel-dominated (3-5%) | Systematic homopolymer length miscalls |
| **PacBio HiFi** | Ultra-low balanced (0.01%) | Minimal position dependence |

---

## Architecture

```
synthdna_lab/
├── config.py              # Company profiles, sequencing backends, pipeline config
├── sequence_generator.py  # 3rd-order Markov chain + synthesis constraint filters
├── secondary_structure.py # Fast hairpin/palindrome detection
├── error_injector.py      # 3-layer error engine (Synthesis → PCR → Sequencing)
├── trace_generator.py     # Negative Binomial coverage model
├── pipeline.py            # Orchestrator + multiprocessing + export (FASTA/CSV/PT)
└── web/
    ├── app.py             # Flask REST API
    ├── templates/         # Jinja2 HTML
    └── static/            # CSS + JS (Canvas charts, DNA animations)
```

### Error Injection Pipeline

```
Original Sequence
       │
       ▼
┌──────────────────┐
│  Layer 1:        │  Position-dependent + homopolymer-biased
│  SYNTHESIS       │  Del:Sub:Ins ratio varies by tech class
│  (vendor-specific)│  Chemical: 60:30:10  |  Enzymatic: 40:40:20
└────────┬─────────┘
         ▼
┌──────────────────┐
│  Layer 2:        │  Very rare (~5e-6/base)
│  PCR             │  KAPA HiFi polymerase model
│  AMPLIFICATION   │  Substitution-dominated
└────────┬─────────┘
         ▼
┌──────────────────┐
│  Layer 3:        │  Backend-specific
│  SEQUENCING      │  Illumina: sub-only + 3' decay
│  (backend-select)│  Nanopore: IDS + homopolymer collapse
└────────┬─────────┘
         ▼
   Corrupted Trace
```

---

## Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/YacineBaghli/SynthDNA-Lab.git
cd SynthDNA-Lab

# Install dependencies
pip install -e .

# Or install with PyTorch export support
pip install -e ".[torch]"
```

### Launch the Web UI

```bash
python -m synthdna_lab.web.app
# Open http://127.0.0.1:5000
```

### Python API

```python
from synthdna_lab import PipelineConfig, TwistRealisticGenerator

# Configure for Twist + Illumina (default)
config = PipelineConfig(target_len=110, dataset_size=10000)
config.set_profile('twist')
config.set_sequencing('illumina')

# Generate a single sample
gen = TwistRealisticGenerator(config)
center, traces = gen.generate_sample(index=0)

print(f"Center: {center[:50]}...")
print(f"Traces: {len(traces)} reads")
for t in traces[:3]:
    print(f"  {t[:50]}... (len={len(t)}, diff={len(center)-len(t)})")
```

### Generate a Full Dataset

```python
from synthdna_lab import PipelineConfig, generate_dataset

config = PipelineConfig(
    target_len=110,
    dataset_size=100_000,
    num_workers=8,
)
config.set_profile('photolitho_ethz')  # ETH photolithographic
config.set_sequencing('nanopore')       # Oxford Nanopore

stats = generate_dataset(
    config=config,
    profile_key='photolitho_ethz',
    output_dir='./datasets',
    output_format='fasta',  # or 'csv', 'pt'
)
print(f"Generated {stats['dataset_size']:,} samples in {stats['generation_time_sec']}s")
```

---

## Running Tests

```bash
pip install -e ".[dev]"
pytest tests/ -v
```

---

## Key Technical Decisions

| Decision | Rationale |
|---|---|
| **3-layer error model** | Preserves distinct physical signatures (synthesis ≠ sequencing). A flat IDS model conflates these. |
| **3rd-order Markov chain** | Human genome trinucleotide frequencies ensure biologically plausible sequence composition. |
| **Nanopore homopolymer collapse** | Models systematic run-length miscalls (not random deletions) — fundamentally different error mechanism. |
| **Negative Binomial coverage** | Captures PCR amplification overdispersion (Poisson is inadequate). |
| **Tech-class-aware substitution matrices** | Chemical synthesis: depurination bias (A↔G). Enzymatic: uniform. |

---

## References

- Grass, R.N., Stark, W.J., Heckel, R. et al. — *Low cost DNA data storage using photolithographic synthesis* (ETH Zurich, 2024)
- Twist Bioscience — Silicon-based synthesis platform specifications
- DNA Script — SYNTAX enzymatic synthesis technical notes
- Ansa Biotechnologies — Tethered-TdT error characterization
- Evonetix — Thermal chip Binary Assembly documentation

---

## License

MIT License — see [LICENSE](LICENSE) for details.

---

<div align="center">

**Built for DNA Data Storage Research**


</div>
