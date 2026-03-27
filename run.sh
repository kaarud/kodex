#!/usr/bin/env bash
# run.sh — Lance le pipeline de filtrage sur un fichier NiourK XLSX
# Usage: bash run.sh <fichier.xlsx> [--beta]
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="ngs_reanalysis"

if [ $# -lt 1 ]; then
    echo "Usage: bash $0 <fichier_niourk.xlsx> [--beta]"
    echo ""
    echo "Commandes disponibles :"
    echo "  bash $0 sample.xlsx                   Filtre complet"
    echo "  bash $0 sample.xlsx --beta            Sans pré-exclusions"
    echo "  bash $0 --aggregate <dossier_xlsx>    Agrège les annotations"
    echo "  bash $0 --test                        Lance les tests"
    exit 1
fi

eval "$(micromamba shell hook -s bash 2>/dev/null || true)"

case "${1:-}" in
    --test)
        echo "==> Tests du pipeline..."
        micromamba run -n "$ENV_NAME" python -m pytest "$SCRIPT_DIR/tests/" -q
        ;;
    --aggregate)
        shift
        TARGET="${1:-.}"
        echo "==> Agrégation des annotations : $TARGET"
        micromamba run -n "$ENV_NAME" python "$SCRIPT_DIR/aggregate_annotations.py" "$TARGET"
        ;;
    *)
        XLSX="$1"
        shift
        echo "==> Filtrage : $XLSX"
        micromamba run -n "$ENV_NAME" python "$SCRIPT_DIR/vcf_filter.py" "$XLSX" --config "$SCRIPT_DIR/config.yaml" "$@"
        ;;
esac
