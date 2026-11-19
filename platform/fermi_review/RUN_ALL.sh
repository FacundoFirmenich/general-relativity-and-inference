#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python repaired/catalog_audit.py
.venv/bin/python repaired/audit_originals.py
.venv/bin/python repaired/audit_hunter.py
printf '%s\n' 'Resultados / Results: results/'
