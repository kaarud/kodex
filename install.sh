#!/usr/bin/env bash
# install.sh — Installation de KODEX sur serveur Ubuntu
# Usage: bash install.sh [--prefix /chemin/micromamba]
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="kodex"
MAMBA_PREFIX="${1:-$HOME/micromamba}"

# ---------------------------------------------------------------
# 1. Installer micromamba si absent
# ---------------------------------------------------------------
if ! command -v micromamba &>/dev/null; then
    echo "==> Installation de micromamba..."
    "${SHELL}" <(curl -L micro.mamba.pm/install.sh) <<EOF
$MAMBA_PREFIX
Y
EOF
    export MAMBA_ROOT_PREFIX="$MAMBA_PREFIX"
    eval "$(micromamba shell hook -s bash)"
    echo "[OK] micromamba installé dans $MAMBA_PREFIX"
else
    echo "[OK] micromamba déjà installé : $(which micromamba)"
    eval "$(micromamba shell hook -s bash)"
fi

# ---------------------------------------------------------------
# 2. Créer l'environnement
# ---------------------------------------------------------------
if micromamba env list | grep -q "$ENV_NAME"; then
    echo "[OK] Environnement '$ENV_NAME' existe déjà — mise à jour..."
    micromamba update -n "$ENV_NAME" -f "$SCRIPT_DIR/environment.yml" -y
else
    echo "==> Création de l'environnement '$ENV_NAME'..."
    micromamba create -f "$SCRIPT_DIR/environment.yml" -y
fi

# ---------------------------------------------------------------
# 3. Vérification
# ---------------------------------------------------------------
echo ""
echo "==> Vérification de l'installation..."
micromamba run -n "$ENV_NAME" python -c "
import sys
print(f'Python  : {sys.version}')
import pandas; print(f'pandas  : {pandas.__version__}')
import openpyxl; print(f'openpyxl: {openpyxl.__version__}')
import yaml; print(f'PyYAML  : {yaml.__version__}')
import requests; print(f'requests: {requests.__version__}')
"

# ---------------------------------------------------------------
# 4. Tests
# ---------------------------------------------------------------
echo ""
echo "==> Exécution des tests..."
cd "$SCRIPT_DIR"
micromamba run -n "$ENV_NAME" python -m pytest tests/ -q

echo ""
echo "======================================="
echo "  Installation terminée avec succès"
echo "  Activation : micromamba activate $ENV_NAME"
echo "======================================="
