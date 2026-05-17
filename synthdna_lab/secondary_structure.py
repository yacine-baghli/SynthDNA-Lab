"""
Secondary Structure Filter — Fast Hairpin Detection.

Uses a palindromic reverse-complement scan to detect potential
hairpin-forming regions that would block reagent access during
silicon-based synthesis.
"""


COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return ''.join(COMPLEMENT.get(b, 'N') for b in reversed(seq))


def has_stable_hairpin(seq: str, min_stem: int = 6,
                       min_loop: int = 3, max_loop: int = 8) -> bool:
    """
    Fast approximate hairpin detection.
    
    Scans for reverse-complement palindromes with:
      - stem length >= min_stem (6 bp → ΔG ≈ -5 kcal/mol)
      - loop length between min_loop and max_loop
    
    Returns True if a stable hairpin is found.
    """
    L = len(seq)
    for i in range(min_stem, L - min_loop - min_stem + 1):
        for loop_len in range(min_loop, min(max_loop + 1, L - i - min_stem + 1)):
            j = i + loop_len
            # Check if stem of length min_stem forms
            stem_left = seq[i - min_stem:i]
            stem_right = seq[j:j + min_stem]
            if reverse_complement(stem_left) == stem_right:
                return True
    return False
