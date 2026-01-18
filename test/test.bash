#!/usr/bin/env bash

set -eo pipefail

PATO_BIN="$1"
shift

MODE="$1"
INDEX="$2"
shift 2

ARGS=()
for arg in "$@"; do
  # Skip empty arguments coming from CMake optional lists
  [[ -z "$arg" ]] && continue

  # Split the argument into words
  for word in $arg; do
    ARGS+=("$word")
  done
done

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
INPUT_DIR="$ROOT_DIR/input"
REF_DIR="$ROOT_DIR/ref"
OUT_DIR="$PWD/output"

mkdir -p "$OUT_DIR"

TFO_FILE="$INPUT_DIR/tfo.fa"
TTS_FILE="$INPUT_DIR/tts.fa"

OUT_PREFIX="$OUT_DIR/${MODE}${INDEX}"

case "$MODE" in
  tfo)
    "$PATO_BIN" "${ARGS[@]}" -ss "$TFO_FILE" -o "$OUT_PREFIX"
    ;;
  tts)
    "$PATO_BIN" "${ARGS[@]}" -ds "$TTS_FILE" -o "$OUT_PREFIX"
    ;;
  tpx)
    "$PATO_BIN" "${ARGS[@]}" -ss "$TFO_FILE" -ds "$TTS_FILE" -o "$OUT_PREFIX"
    ;;
  *)
    echo "Unknown test mode: $MODE" >&2
    exit 1
    ;;
esac

sort "${OUT_PREFIX}.out" > "${OUT_PREFIX}.sorted.out"
sort "${OUT_PREFIX}.summary" > "${OUT_PREFIX}.sorted.summary"
sort "$REF_DIR/${MODE}${INDEX}.out" > "$OUT_DIR/${MODE}${INDEX}.ref.sorted.out"
sort "$REF_DIR/${MODE}${INDEX}.summary" > "$OUT_DIR/${MODE}${INDEX}.ref.sorted.summary"

diff -u "$OUT_DIR/${MODE}${INDEX}.ref.sorted.out" "${OUT_PREFIX}.sorted.out"
diff -u "$OUT_DIR/${MODE}${INDEX}.ref.sorted.summary" "${OUT_PREFIX}.sorted.summary"

exit 0
