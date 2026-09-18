#!/usr/bin/env bash
# Build and run the VALIDATION-ONLY closing sweep of the XC table generator
# against the reference tables in data/tables/fortran_native.
#
#   scripts/sweep_reference_tables.sh data/tables/fortran_native/xc_table_u*.dat
#
# Delete this script together with original/ and data/tables/.
#
# The sources are compiled here instead of linking fpm's archive on purpose:
# fpm 0.12.0 cannot carry the OpenMP flags in fpm.toml (see the OpenMP note
# there), so the build directory holds several profiles at once and picking one
# by timestamp silently linked the *debug* library, which made the reported
# generation times a factor of six too large.  The file list below is the
# module dependency order of `sweep_reference_tables.f90`; extend it if that
# program grows a new `use`.
#
# Kept out of fpm.toml deliberately: the sweep takes tens of minutes and must
# never run as part of `fpm test`.
set -euo pipefail

cd "$(dirname "$0")/.."

OUT=${TMPDIR:-/tmp}/lsdaks_sweep
mkdir -p "$OUT/mod"

gfortran -O3 -march=native -fopenmp -J "$OUT/mod" \
    src/types/lsda_constants.f90 \
    src/types/lsda_errors.f90 \
    src/bethe_ansatz/table_io.f90 \
    src/bethe_ansatz/lieb_wu_integral.f90 \
    src/bethe_ansatz/bethe_tables.f90 \
    src/xc_functional/spline2d.f90 \
    src/xc_functional/xc_lsda.f90 \
    scripts/sweep_reference_tables.f90 \
    -llapack -lblas -o "$OUT/sweep"

exec "$OUT/sweep" "$OUT" "$@"
