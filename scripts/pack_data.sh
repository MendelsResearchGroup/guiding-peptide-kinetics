#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
OUT_DIR="${1:-$REPO_ROOT/data_archives}"

if ! command -v zip >/dev/null 2>&1; then
  echo "'zip' is required but not found on PATH." >&2
  exit 1
fi

mkdir -p "$OUT_DIR"

# Rebuild canonical archive names so stale content does not persist.
rm -f "$OUT_DIR"/data_core.zip
rm -f "$OUT_DIR"/hlda_trajectories_*.zip

pushd "$REPO_ROOT" >/dev/null

# Core data bundle: everything in data/ except large trajectory cache.
zip -r -q "$OUT_DIR/data_core.zip" data -x 'data/hlda_trajectories/*'

# One archive per trajectory subdirectory to keep files manageable.
if [[ -d data/hlda_trajectories ]]; then
  shopt -s nullglob
  for d in data/hlda_trajectories/*; do
    [[ -d "$d" ]] || continue
    name="$(basename "$d")"
    zip -r -q "$OUT_DIR/hlda_trajectories_${name}.zip" "$d"
  done
fi

popd >/dev/null

count=$(find "$OUT_DIR" -maxdepth 1 -type f -name '*.zip' | wc -l)
echo "Created $count archive(s) in $OUT_DIR"
