#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

ARCHIVE_DIR="${1:-}"
if [[ -z "$ARCHIVE_DIR" ]]; then
  if [[ -d "$REPO_ROOT/data_archives" ]]; then
    ARCHIVE_DIR="$REPO_ROOT/data_archives"
  elif [[ -d "$REPO_ROOT/data/archives" ]]; then
    ARCHIVE_DIR="$REPO_ROOT/data/archives"
  else
    echo "No archive directory found. Expected one of:" >&2
    echo "  $REPO_ROOT/data_archives" >&2
    echo "  $REPO_ROOT/data/archives" >&2
    echo "Or pass a directory explicitly: ./scripts/unpack_data.sh <archive_dir>" >&2
    exit 1
  fi
fi

if ! command -v unzip >/dev/null 2>&1; then
  echo "'unzip' is required but not found on PATH." >&2
  exit 1
fi

shopt -s nullglob
zips=("$ARCHIVE_DIR"/*.zip)
if (( ${#zips[@]} == 0 )); then
  echo "No zip archives found in $ARCHIVE_DIR" >&2
  exit 1
fi

echo "Unpacking ${#zips[@]} archive(s) from $ARCHIVE_DIR into $REPO_ROOT"
for z in "${zips[@]}"; do
echo "  -> $(basename "$z")"
  unzip -n "$z" -d "$REPO_ROOT"
done
