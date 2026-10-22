"""Neighbour-gene windows must stay in GRCh38.

Both entries of `_NEIGHBOR_GENES` once carried hg19 remnants under a comment
announcing GRCh38: PKD1 kept its hg19 end coordinate (2 185 899) while its
start had been converted, and RHEB sat ~32 kb off. On RHEB, only one of the
eight amplicons the panel actually sequences fell inside the declared window.

The bug was silent — it only affects the blue "neighbour gene" marking, never
the filtering — which is precisely why it survived. These tests freeze the
Ensembl GRCh38 gene bounds so a future edit has to be deliberate.
"""
from vcf_filter import _NEIGHBOR_GENES, _in_neighbor_gene

# Ensembl REST, GRCh38, whole-gene bounds (1-based, no margin)
EXPECTED = {
    "PKD1": ("chr16", 2088708, 2135898),
    "RHEB": ("chr7", 151465976, 151520120),
}

# hg19 bounds — must never reappear in the GRCh38 table
HG19_BOUNDS = {
    "PKD1": ("chr16", 2138711, 2185899),
    "RHEB": ("chr7", 151163098, 151217206),
}


def test_bounds_are_grch38():
    for gene, (chrom, start, end) in EXPECTED.items():
        entry = _NEIGHBOR_GENES[gene]
        assert (entry["chr"], entry["start"], entry["end"]) == (chrom, start, end)


def test_no_hg19_remnant():
    """The old defect was a *mixed* window — one converted bound, one not."""
    for gene, (_, hg19_start, hg19_end) in HG19_BOUNDS.items():
        entry = _NEIGHBOR_GENES[gene]
        assert entry["start"] != hg19_start, f"{gene}: start is the hg19 bound"
        assert entry["end"] != hg19_end, f"{gene}: end is the hg19 bound"


def test_rheb_covers_the_sequenced_gene():
    """The whole RHEB gene body must fall inside the window, not one amplicon."""
    for pos in (151465976, 151490000, 151520120):
        assert _in_neighbor_gene("chr7", pos) == "RHEB"


def test_pkd1_window_stops_at_the_gene_end():
    assert _in_neighbor_gene("chr16", 2135898) == "PKD1"
    # 2 185 899 is the hg19 end — in GRCh38 it lies well past PKD1
    assert _in_neighbor_gene("chr16", 2185899) is None


def test_outside_returns_none():
    assert _in_neighbor_gene("chr16", 1000000) is None
    assert _in_neighbor_gene("chr1", 2100000) is None
