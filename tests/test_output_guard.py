"""A *_filtered.xlsx holds hand-typed USER_CLASS/ANNOT/COM: never overwrite it.

vcf_filter opened the output in write mode, so re-running it on a sample wiped
weeks of review with no warning and no backup.
"""

import pytest
from vcf_filter import resolve_output_path


def test_plain_name_when_nothing_exists(tmp_path):
    src = tmp_path / "run_2_Singularity_100007_S6.xlsx"
    src.touch()
    assert resolve_output_path(src).name == "run_2_Singularity_100007_S6_filtered.xlsx"


def test_existing_output_is_never_clobbered(tmp_path):
    src = tmp_path / "run_2_Singularity_100007_S6.xlsx"
    src.touch()
    (tmp_path / "run_2_Singularity_100007_S6_filtered.xlsx").touch()
    out = resolve_output_path(src)
    assert out.name == "run_2_Singularity_100007_S6_filtered_2.xlsx"
    assert not out.exists()


def test_suffix_keeps_climbing(tmp_path):
    src = tmp_path / "S6.xlsx"
    src.touch()
    for name in ("S6_filtered.xlsx", "S6_filtered_2.xlsx", "S6_filtered_3.xlsx"):
        (tmp_path / name).touch()
    assert resolve_output_path(src).name == "S6_filtered_4.xlsx"


def test_force_reclaims_the_plain_name(tmp_path):
    src = tmp_path / "S6.xlsx"
    src.touch()
    (tmp_path / "S6_filtered.xlsx").touch()
    assert resolve_output_path(src, force=True).name == "S6_filtered.xlsx"


def test_incremented_name_stays_out_of_the_aggregation_glob(tmp_path):
    """'*_filtered.xlsx' must not match the new file: the annotated original
    stays the one aggregate_annotations reads, and no sample is counted twice."""
    src = tmp_path / "S6.xlsx"
    src.touch()
    (tmp_path / "S6_filtered.xlsx").touch()
    out = resolve_output_path(src)
    assert out not in tmp_path.glob("*_filtered.xlsx")
