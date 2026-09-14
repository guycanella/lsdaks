#!/usr/bin/env bash
# Compila o código C++ original (original/) fora da árvore, em build/cpp/,
# e prepara o diretório para rodar casos de referência no macOS.
#
# Uso:
#   scripts/build_cpp_reference.sh            # compila e prepara build/cpp
#   scripts/build_cpp_reference.sh --run      # idem e roda o caso de fumaça (uniforme, U=4)
#
# Por que não usar o makefile original: ele usa -static e -pg, que o clang do
# macOS não suporta, e gera os .o dentro de original/, que é referência imutável.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRC="$ROOT/original"
OUT="$ROOT/build/cpp"
BIN="$OUT/lsdaks_cpp"

CXX="${CXX:-clang++}"
CXXFLAGS="${CXXFLAGS:--O2 -std=c++11 -Wno-deprecated}"

mkdir -p "$OUT"

echo "==> Compilando com $CXX $CXXFLAGS"
"$CXX" $CXXFLAGS -o "$BIN" \
  "$SRC/lsdaks.cc" "$SRC/lsda_potential.cc" "$SRC/spline2D.cc" "$SRC/giv.cc" \
  "$SRC/linear.cc" "$SRC/lsda_stop.cc" "$SRC/lsda_simetria.cc" "$SRC/lsda_twist.cc" \
  "$SRC/lsda_interface.cc" -lm
echo "    binário: $BIN"

echo "==> Linkando tabelas e potential_param.dat em $OUT"
ln -sf "$SRC/potential_param.dat" "$OUT/potential_param.dat"
for f in "$SRC"/lsda_hub_u*; do
  ln -sf "$f" "$OUT/$(basename "$f")"
done

# O programa lê 11 linhas do stdin: Na Nup Ndn bc u pot symm funcional filename interative [phase] TOL?
# Ver original/lsda_interface.cc:78-161. Se bc=2, a fase (em unidades de pi) vem antes da resposta do TOL.
cat > "$OUT/input_uniform_u4.dat" <<'EOF'
90
45
45
0
4
1
0
n
ref_uniform_u4
0
y
1e-10
EOF

if [[ "${1:-}" == "--run" ]]; then
  echo "==> Rodando caso de fumaça: L=90, N=45/45, OBC, U=4, potencial uniforme"
  (cd "$OUT" && time ./lsdaks_cpp < input_uniform_u4.dat > run_uniform_u4.log 2>&1)
  echo "==> Resultado:"
  sed -n '1,13p' "$OUT/ref_uniform_u4"
fi

echo
echo "Pronto. Para rodar outro caso:"
echo "  cd $OUT && ./lsdaks_cpp < <seu_input>.dat"
echo "Saída: arquivo com o nome dado na linha 9 do input, mais evolution.dat (iter, E_total, Mix)."
