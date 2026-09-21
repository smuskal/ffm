#!/usr/bin/env bash
# Thin wrapper so the tool works without installing anything into the system.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
[ -d "$HERE/.venv" ] && . "$HERE/.venv/bin/activate"
exec python3 -m familyfm.cli "$@"
