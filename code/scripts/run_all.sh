#!/usr/bin/env bash
set -euo pipefail
OUT=out; mkdir -p "$OUT"
threads="${1:-8}"
for r in Rule/*.txt; do
  b="$(basename "$r" .txt)"
  lay="case/${b}.txt"
  out="${OUT}/${b}.txt"
  echo "==> $b"
  ./trace -layout "$lay" -rule "$r" -output "$out" -thread "$threads"
  if [ -x ./AnswerChecker ]; then
    ./AnswerChecker -std "answer/${b}.txt" -out "$out"
  fi
done
