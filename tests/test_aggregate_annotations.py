"""Tests for aggregate_annotations._collect_annotated row selection."""

from aggregate_annotations import _collect_annotated


def _rows(records):
    """Wrap records with distinct CHR/POS/REF/ALT so dedup keeps them separate."""
    out = []
    for i, rec in enumerate(records):
        out.append({"CHR": "chr9", "POS": 100 + i, "REF": "A", "ALT": "T", **rec})
    return out


def test_collect_keeps_sure_and_doubt_scores(make_filtered_xlsx):
    """Plain numbers, a dot, and any value ending in '?' (doubt) are extracted."""
    path = make_filtered_xlsx({
        "HIGH_IMPACT": _rows([
            {"USER_CLASS": "3"},     # sure
            {"USER_CLASS": "."},     # sure
            {"USER_CLASS": "0?"},    # doubt — maybe artefact, keep for review
            {"USER_CLASS": "path?"}, # doubt (word + ?)
        ]),
    })
    kept = _collect_annotated(path)
    assert {r["USER_CLASS"] for r in kept} == {"3", ".", "0?", "path?"}


def test_collect_drops_bare_zero_artefact(make_filtered_xlsx):
    """A bare USER_CLASS '0' (confirmed artefact) with no annotation is dropped."""
    path = make_filtered_xlsx({
        "HIGH_IMPACT": _rows([
            {"USER_CLASS": "0"},          # artefact, no note → drop
            {"USER_CLASS": " 0 "},        # whitespace around → still artefact → drop
            {"USER_CLASS": "5"},          # keep
        ]),
    })
    kept = _collect_annotated(path)
    assert [str(r["USER_CLASS"]) for r in kept] == ["5"]


def test_collect_keeps_zero_artefact_when_annotated(make_filtered_xlsx):
    """A '0' artefact still surfaces if the reviewer wrote a USER_ANNOT note."""
    path = make_filtered_xlsx({
        "HIGH_IMPACT": _rows([
            {"USER_CLASS": "0", "USER_ANNOT": "strand bias, ignore"},
            {"USER_CLASS": "0", "USER_ANNOT": ""},
        ]),
    })
    kept = _collect_annotated(path)
    assert len(kept) == 1
    assert kept[0]["USER_ANNOT"] == "strand bias, ignore"


# --- _collect_all: every variant, scored, for sensitivity control -----------

def test_collect_all_keeps_unannotated_rows(make_filtered_xlsx):
    """Unlike _collect_annotated, --all keeps rows with no USER_* at all."""
    from aggregate_annotations import _collect_all
    path = make_filtered_xlsx({
        "HIGH_IMPACT": _rows([
            {"SYMBOL": "TSC1", "USER_CLASS": "3"},
            {"SYMBOL": "TSC2"},                      # never reviewed
            {"SYMBOL": "TSC2", "USER_CLASS": "0"},   # confirmed artefact
        ]),
    })
    assert len(_collect_all(path)) == 3


def test_collect_all_marks_which_rows_were_retained(make_filtered_xlsx):
    """RETENU flags the rows _collect_annotated would have kept, and only those."""
    from aggregate_annotations import _collect_all
    path = make_filtered_xlsx({
        "HIGH_IMPACT": _rows([
            {"SYMBOL": "TSC1", "USER_CLASS": "3"},    # retained
            {"SYMBOL": "TSC1", "USER_ANNOT": "note"}, # retained
            {"SYMBOL": "TSC2"},                       # not reviewed
            {"SYMBOL": "TSC2", "USER_CLASS": "0"},    # bare 0 = artefact
        ]),
    })
    flags = [r["RETENU"] for r in _collect_all(path)]
    assert flags == ["oui", "oui", "", ""]


def test_collect_all_excludes_extra_genes(make_filtered_xlsx):
    """RHEB/PKD1 are kept in the workbook but stay out of the ranking."""
    from aggregate_annotations import _collect_all
    path = make_filtered_xlsx({
        "HIGH_IMPACT": _rows([{"SYMBOL": "TSC1"}, {"SYMBOL": "TSC2"}]),
        "RHEB_PKD1": _rows([{"SYMBOL": "RHEB"}, {"SYMBOL": "PKD1"}]),
    })
    kept = _collect_all(path, extra_genes=["RHEB", "PKD1"])
    assert {r["SYMBOL"] for r in kept} == {"TSC1", "TSC2"}


def test_collect_all_dedups_across_sheets(make_filtered_xlsx):
    """One variant in two sheets is one row, both sheet names recorded."""
    from aggregate_annotations import _collect_all
    row = {"CHR": "chr9", "POS": 500, "REF": "A", "ALT": "T", "SYMBOL": "TSC1"}
    path = make_filtered_xlsx({"HIGH_IMPACT": [row], "FQ_LOW": [row]})
    kept = _collect_all(path)
    assert len(kept) == 1
    assert set(kept[0]["SOURCE_SHEET"].split(" | ")) == {"HIGH_IMPACT", "FQ_LOW"}


def test_sort_by_score_can_rank_across_patients():
    """--all ranks globally: the best variant leads whoever it belongs to."""
    import pandas as pd
    from aggregate_annotations import _sort_by_score
    df = pd.DataFrame([
        {"SAMPLE": "S1", "SCORE": 0.2},
        {"SAMPLE": "S2", "SCORE": 0.9},
        {"SAMPLE": "S1", "SCORE": 0.5},
    ])
    ranked = _sort_by_score(df, group_by_sample=False)
    assert list(ranked["SCORE"]) == [0.9, 0.5, 0.2]


def test_aggregate_all_ranks_on_a_real_score(make_filtered_xlsx):
    """--all is worthless without SCORE: ranking must not fall back to position."""
    from aggregate_annotations import aggregate_all
    path = make_filtered_xlsx({
        "HIGH_IMPACT": [
            # Weak first, strong second: a positional sort would keep this order.
            {"CHR": "chr9", "POS": 100, "REF": "A", "ALT": "T", "SYMBOL": "TSC1",
             "IMPACT": "MODIFIER", "CADD_PHRED": 1.0, "CLINVAR": None},
            {"CHR": "chr9", "POS": 200, "REF": "A", "ALT": "T", "SYMBOL": "TSC1",
             "IMPACT": "HIGH", "CADD_PHRED": 35.0, "CLINVAR": "Pathogenic"},
        ],
    })
    df = aggregate_all(path.parent, files=[path], every_variant=True)
    assert "SCORE" in df.columns
    assert list(df["POS"]) == [200, 100]


def test_partial_score_column_is_recomputed(make_filtered_xlsx, tmp_path):
    """Workbooks from different KODEX versions mix: some carry SCORE, some not.

    Concatenating them yields a SCORE column that exists but is mostly NaN, and
    testing `"SCORE" not in df.columns` skips the recompute — the unscored rows
    then sort last, burying exactly the variants under review.
    """
    import pandas as pd
    from aggregate_annotations import aggregate_all

    def book(name, with_score):
        p = tmp_path / name
        rows = [{"CHR": "chr9", "POS": 100, "REF": "A", "ALT": "T", "SYMBOL": "TSC1",
                 "IMPACT": "HIGH", "CADD_PHRED": 35.0, "DP": 500, "AF_PCT": 45.0,
                 "SB": 0.5, "CALLNB": 4, "RNCNT_TOTAL": 0, "PJCNT_TOTAL": 0}]
        if with_score:
            rows[0]["SCORE"] = 0.42
        with pd.ExcelWriter(p, engine="openpyxl") as w:
            pd.DataFrame({"Key": ["date"], "Value": ["x"]}).to_excel(
                w, sheet_name="_summary", index=False)
            pd.DataFrame(rows).to_excel(w, sheet_name="HIGH_IMPACT", index=False)
        return p

    files = [book("run_1_111_S1_filtered.xlsx", True),
             book("run_1_222_S2_filtered.xlsx", False)]
    df = aggregate_all(tmp_path, files=files, every_variant=True)
    assert df["SCORE"].notna().all(), "toute ligne doit porter un score"
