#!/bin/bash
# Reproducible verification for the UID-3619 decoder hardening.
# Builds an ASan+UBSan static library and a standalone harness, then runs every
# regression seed and every malformed PoC and reports crash / sanitizer counts.
set -u

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="$ROOT/build-verify"

echo "==> configuring ASan+UBSan static build"
cmake -S "$ROOT" -B "$BUILD" -DCMAKE_BUILD_TYPE=Debug -DBUILD_SHARED_LIBS=OFF \
  -DWAVPACK_BUILD_PROGRAMS=OFF -DWAVPACK_ENABLE_LEGACY=ON \
  -DWAVPACK_INSTALL_DOCS=OFF -DWAVPACK_INSTALL_CMAKE_MODULE=OFF \
  -DWAVPACK_INSTALL_PKGCONFIG_MODULE=OFF \
  -DCMAKE_C_FLAGS="-g3 -O1 -fsanitize=address,undefined -fno-sanitize=signed-integer-overflow -fno-omit-frame-pointer"
cmake --build "$BUILD" -j"$(nproc)"

echo "==> building standalone harness"
gcc -g3 -O1 -fsanitize=address,undefined -fno-sanitize=signed-integer-overflow \
  -fno-omit-frame-pointer -I"$ROOT/include" "$ROOT/fuzz-local/harness.c" \
  "$BUILD/libwavpack.a" -lm -lpthread -o "$ROOT/fuzz-local/harness-verify"

export ASAN_OPTIONS=detect_leaks=0:abort_on_error=1
export UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1

total=0; crashed=0
echo "==> running all regression seeds and malformed PoCs"
for f in "$ROOT"/fuzzing/regression/*.wv "$ROOT"/artifacts/malformed-samples/*.wv; do
  total=$((total + 1))
  if ! "$ROOT/fuzz-local/harness-verify" "$f" >/dev/null 2>"$ROOT/fuzz-local/.last.err"; then
    crashed=$((crashed + 1))
    echo "CRASH: $f"
    head -3 "$ROOT/fuzz-local/.last.err"
  fi
done

echo "----------------------------------------"
echo "files tested : $total"
echo "crashes      : $crashed"
echo "sanitizer    : ASan + UBSan (signed integer overflow excluded, matching upstream)"
if [ "$crashed" -eq 0 ]; then
  echo "RESULT       : ZERO CRASHES"
  exit 0
fi
exit 1
