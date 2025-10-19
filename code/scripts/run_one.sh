#!/usr/bin/env bash
set -euo pipefail
if [ $# -lt 3 ]; then
  echo "Usage: $0 <layout.txt> <rule.txt> <out.txt> [threads]"
  exit 1
fi
threads="${4:-8}"
./trace -layout "$1" -rule "$2" -output "$3" -thread "$threads"
