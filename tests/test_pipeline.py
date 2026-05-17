"""Basic smoke tests for the SynthDNA Lab pipeline."""

import pytest
from synthdna_lab.config import PipelineConfig, COMPANY_PROFILES, SEQUENCING_BACKENDS
from synthdna_lab.pipeline import TwistRealisticGenerator


@pytest.fixture(params=list(COMPANY_PROFILES.keys()))
def profile_key(request):
    return request.param


@pytest.fixture(params=list(SEQUENCING_BACKENDS.keys()))
def seq_backend(request):
    return request.param


class TestSequenceGeneration:
    """Test that all profiles generate valid sequences."""

    def test_generates_center_sequence(self, profile_key):
        config = PipelineConfig(target_len=110)
        config.set_profile(profile_key)
        gen = TwistRealisticGenerator(config)
        center, traces = gen.generate_sample(0, seed_offset=42)
        assert len(center) > 0
        assert all(b in 'ACGT' for b in center)

    def test_generates_traces(self, profile_key):
        config = PipelineConfig(target_len=110)
        config.set_profile(profile_key)
        gen = TwistRealisticGenerator(config)
        center, traces = gen.generate_sample(0, seed_offset=42, min_traces=3, max_traces=10)
        assert len(traces) >= 3
        for t in traces:
            assert len(t) > 0
            assert all(b in 'ACGT' for b in t)

    def test_gc_within_range(self, profile_key):
        config = PipelineConfig(target_len=110)
        config.set_profile(profile_key)
        gen = TwistRealisticGenerator(config)
        for i in range(20):
            center, _ = gen.generate_sample(i, seed_offset=42)
            gc = (center.count('G') + center.count('C')) / len(center)
            low, high = config.profile.gc_range
            assert low <= gc <= high, f"GC {gc:.3f} outside [{low}, {high}]"


class TestSequencingBackends:
    """Test that sequencing backends don't crash."""

    def test_backend_produces_traces(self, seq_backend):
        config = PipelineConfig(target_len=110)
        config.set_profile('twist')
        config.set_sequencing(seq_backend)
        gen = TwistRealisticGenerator(config)
        center, traces = gen.generate_sample(0, seed_offset=42)
        assert len(traces) >= 1


class TestPhotolithographic:
    """Photolithographic profile should show dramatically higher error rates."""

    def test_high_error_rate(self):
        config = PipelineConfig(target_len=110)
        config.set_profile('photolitho_ethz')
        gen = TwistRealisticGenerator(config)

        total_len_diff = 0
        n = 50
        for i in range(n):
            center, traces = gen.generate_sample(i, seed_offset=42)
            for t in traces:
                total_len_diff += abs(len(center) - len(t))

        # Photolithographic should have much larger average length diffs
        avg_diff = total_len_diff / n
        assert avg_diff > 5, f"Expected large length diffs for photolitho, got {avg_diff}"
