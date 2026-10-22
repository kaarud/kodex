"""Le filtre doit retrouver le variant qui a motive la mesure du bruit de fond."""

import bg_filter as B

# Extrait reel du commentaire ecrit par run_background.py pour
# TSC2 c.4090dup chez 100015_S14 (run 1).
TEMOIN = """chr16:2084311 C>CT — 100015_S14 (patient)
  DP=2137  A:0 C:2135 G:0 T:2  DEL:0 INS:55
  -> alt 2.57 %

Run — 22 autres echantillons
  100003_S7          DP=1515  A:0 C:1515 G:0 T:0  DEL:0 INS:0  -> 0.00 %

  mediane 0.00 % | max 0.00 % | avec alt : 0/22
  contraste x46.8
  AF hors supplementaires : 2.59 %  (rapport 1.01 — reads primaires — OK)"""


def test_parse_temoin():
    m = B.parse_comment(TEMOIN)
    assert m["AF_BAM_PCT"] == 2.57
    assert m["DP_BAM"] == 2137          # celui du patient, pas du premier autre
    assert m["BG_MEDIAN_PCT"] == 0.0
    assert (m["BG_N_WITH_ALT"], m["BG_N_MEASURED"]) == (0, 22)
    assert m["BG_CONTRAST"] == 46.8
    assert m["BG_PRIMARY_RATIO"] == 1.01


def test_temoin_retenu():
    assert B.keeps(B.parse_comment(TEMOIN))


def test_mesure_absente_nest_pas_bruit_nul():
    """Un site non couvert ne doit pas passer pour un site propre."""
    texte = TEMOIN.replace("mediane 0.00 %", "mediane — %") \
                  .replace("contraste x46.8", "contraste non mesure")
    m = B.parse_comment(texte)
    assert m["BG_MEDIAN_PCT"] is None and m["BG_CONTRAST"] is None
    assert not B.keeps(m)


def test_artefact_de_structure_rejete():
    m = B.parse_comment(TEMOIN.replace("rapport 1.01", "rapport 0.12"))
    assert not B.keeps(m)


def test_bruit_partage_rejete():
    m = B.parse_comment(TEMOIN.replace("mediane 0.00 %", "mediane 0.40 %")
                              .replace("avec alt : 0/22", "avec alt : 19/22")
                              .replace("contraste x46.8", "contraste x1.2"))
    assert not B.keeps(m)


def test_profondeur_faible_rejetee():
    """Le plancher 1/profondeur gonfle le contraste sur un site peu couvert."""
    m = B.parse_comment(TEMOIN.replace("DP=2137", "DP=40"))
    assert not B.keeps(m)
