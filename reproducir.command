#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
exec bash reproducir.bash "${1:-all}"
