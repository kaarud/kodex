"""minibam.py — logique pure : régions, nommage, commande samtools, manifest."""

from pathlib import Path

import pandas as pd
import pytest

from minibam import (
    build_view_cmd,
    merge_regions,
    read_annotations,
    slice_name,
    write_manifest,
)


def test_regions_are_flanked_and_prefixed():
    assert merge_regions([("9", 10_000)], flank=500) == ["chr9:9500-10500"]


def test_regions_clamp_start_to_one():
    assert merge_regions([("chr9", 100)], flank=500) == ["chr9:1-600"]


def test_overlapping_variants_share_one_region():
    # 200 pb d'écart avec un flanc de 500 → une seule tranche
    assert merge_regions([("chr9", 10_000), ("chr9", 10_200)], flank=500) == [
        "chr9:9500-10700"
    ]


def test_distant_variants_stay_separate():
    regions = merge_regions([("chr9", 10_000), ("chr9", 50_000)], flank=500)
    assert regions == ["chr9:9500-10500", "chr9:49500-50500"]


def test_duplicate_positions_collapse():
    assert merge_regions([("chr9", 10_000), ("chr9", 10_000)], flank=500) == [
        "chr9:9500-10500"
    ]


def test_regions_grouped_by_chromosome():
    regions = merge_regions([("chr16", 2_000), ("chr9", 10_000)], flank=100)
    assert regions == ["chr16:1900-2100", "chr9:9900-10100"]


def test_slice_name_keeps_track_type():
    assert slice_name(Path("/x/100003_S7.bam")) == "100003_S7.mini.bam"
    assert slice_name(Path("/x/100003_S7.chim.bam")) == "100003_S7.chim.mini.bam"
    assert slice_name(Path("/x/100003_S7.mutect2.bam")) == "100003_S7.mutect2.mini.bam"


def test_slice_name_flags_the_comparison_track():
    assert slice_name(Path("/x/100010_S1.bam"), "REF_S1") == "REF_100010_S1.mini.bam"


def test_view_cmd_passes_every_region():
    cmd = build_view_cmd(Path("/x/s.bam"), ["chr9:1-2", "chr16:3-4"], Path("/o/s.mini.bam"))
    assert cmd == ["samtools", "view", "-b", "-h", "-o", "/o/s.mini.bam",
                   "/x/s.bam", "chr9:1-2", "chr16:3-4"]


def test_manifest_is_tsv_with_stable_columns(tmp_path):
    entry = {"sample": "s1", "track": "RAW", "file": "s1/s1.mini.bam",
             "source_bam": "/x/s1.bam", "n_variants": 3, "n_regions": 2,
             "bytes": 1234}
    out = tmp_path / "manifest.tsv"
    write_manifest([entry], out)
    df = pd.read_csv(out, sep="\t")
    assert list(df.columns) == ["sample", "track", "file", "source_bam",
                               "n_variants", "n_regions", "bytes"]
    assert df["bytes"].iloc[0] == 1234


def _write_xlsx(path, rows, sheet="ANNOTATIONS"):
    pd.DataFrame(rows).to_excel(path, sheet_name=sheet, index=False)
    return path


def test_read_annotations_rejects_wrong_workbook(tmp_path):
    path = _write_xlsx(tmp_path / "x.xlsx", [{"SAMPLE": "s1"}], sheet="AUTRE")
    with pytest.raises(SystemExit, match="ANNOTATIONS"):
        read_annotations(path)


def test_read_annotations_rejects_missing_columns(tmp_path):
    path = _write_xlsx(tmp_path / "x.xlsx", [{"SAMPLE": "s1", "CHR": "chr9"}])
    with pytest.raises(SystemExit, match="POS"):
        read_annotations(path)


def test_read_annotations_drops_incomplete_rows(tmp_path):
    path = _write_xlsx(tmp_path / "x.xlsx", [
        {"SAMPLE": "s1", "CHR": "chr9", "POS": 10},
        {"SAMPLE": "s1", "CHR": "chr9", "POS": None},
    ])
    assert len(read_annotations(path)) == 1
