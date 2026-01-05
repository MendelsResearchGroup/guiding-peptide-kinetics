#!/bin/bash
set -e

# Usage: ./migrate_outputs_shortnames.sh
# Moves data/<LONG>/output -> data/<SHORT>/output based on src/common/consts.py.
# If data/<SHORT>/output already exists, long outputs are moved into
# data/<SHORT>/output/from_<LONG>_<timestamp> to avoid overwrites.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

cd "$REPO_ROOT"

mapfile -t mappings < <(
  python - <<'PY'
import importlib.util
from pathlib import Path

root = Path('.').resolve()
consts_path = root / 'src' / 'common' / 'consts.py'
spec = importlib.util.spec_from_file_location('consts', consts_path)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
for long, short in mod.long_to_short.items():
    print(f"{long} {short}")
PY
)

timestamp=$(date +%Y%m%d_%H%M%S)

for line in "${mappings[@]}"; do
  long=${line%% *}
  short=${line##* }

  long_dir="data/${long}"
  long_output="${long_dir}/output"
  short_dir="data/${short}"
  short_output="${short_dir}/output"

  if [ ! -d "$long_output" ]; then
    continue
  fi

  if [ ! -d "$short_dir" ]; then
    mkdir -p "$short_dir"
  fi

  if [ ! -d "$short_output" ]; then
    echo "Moving ${long_output} -> ${short_output}"
    mv "$long_output" "$short_output"
    rmdir "$long_dir" 2>/dev/null || true
    continue
  fi

  dest="${short_output}/from_${long}_${timestamp}"
  echo "Moving ${long_output} -> ${dest}"
  mv "$long_output" "$dest"
  rmdir "$long_dir" 2>/dev/null || true
done
