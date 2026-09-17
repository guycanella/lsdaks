#!/usr/bin/env bash
#
# run_examples.sh - smoke test of the input layer of lsdaks.
#
# Every file in examples/ must be accepted by the parser and reach the SCF
# cycle; a malformed input must be REJECTED with a message. Both halves are
# checked here, because the project has no process-level test harness and these
# failures only show up when the executable is actually run:
#
#   0. each examples/*.txt is run EXACTLY AS IT IS IN THE REPOSITORY, byte for
#      byte, and must reach the SCF cycle. This pass exists because an earlier
#      version of this script only ran the rewritten copies produced by
#      prepare_input below, and `awk` silently appends a final newline: that
#      hid a parser bug which rejected every input file whose last line was not
#      newline-terminated (4 of the 6 examples, and input.txt). Any rewriting
#      of the input can mask a byte-level defect, so the originals are tested
#      first and without any preprocessing;
#   1. each examples/*.txt runs again from a rewritten copy, with max_iter
#      forced to 5 so the whole suite takes a second, and gets as far as
#      "Starting Kohn-Sham SCF Cycle";
#   2. an unknown namelist key is rejected with a message naming the key;
#   3. the key `distribution`, removed from &potential, is rejected with the
#      dedicated migration hint instead of a bare "cannot match" message;
#   4. a malformed value (U = 1.2.3) is rejected;
#   5. a genuinely ABSENT namelist group is accepted, with a note (the four
#      groups are optional; only a broken one is an error);
#   6. a run that converges but cannot write its output files exits non-zero
#      (a valid calculation whose record was lost must not report success);
#   7. generate_xc_table never leaves an invalid table behind: either it exits
#      0 and the file it wrote has no NaN/Inf, or it exits non-zero and wrote
#      nothing. Phrased as an invariant rather than "must fail" so it stays
#      valid when the generator is completed (phase 4.5, T21); today the
#      finite-L solver does not converge on part of the U=4 grid and the
#      executable must take the second branch.
#
# The executable is launched from a throwaway directory, including the pass
# that consumes the examples byte-for-byte.  Thus their committed output_prefix
# values cannot overwrite an existing result in the checkout. Exit status: 0
# if every check passed, 1 otherwise.
#
# Usage: scripts/run_examples.sh

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

cd "${PROJECT_DIR}" || exit 1

WORK_DIR="$(mktemp -d /tmp/lsdaks_examples.XXXXXX)"
trap 'rm -rf "${WORK_DIR}"' EXIT

# fpm itself needs the project directory, but the program must not inherit it
# as its working directory: the examples carry ordinary relative output
# prefixes.  Build once, select the executable just produced, then run it from
# WORK_DIR with an absolute table directory.
if ! fpm build >/dev/null; then
    printf 'Cannot build lsdaks.\n' >&2
    exit 1
fi
# Ask fpm for the executable selected by its *current* debug/profile settings.
# Selecting the newest build/*/app/lsdaks is wrong when a release build happens
# to be newer than the default debug target.
EXECUTABLE_REL="$(fpm run lsdaks --runner echo | tail -n 1)"
case "${EXECUTABLE_REL}" in
    /*) EXECUTABLE="${EXECUTABLE_REL}" ;;
    *)  EXECUTABLE="${PROJECT_DIR}/${EXECUTABLE_REL}" ;;
esac
if [ -z "${EXECUTABLE}" ] || [ ! -x "${EXECUTABLE}" ]; then
    printf 'Cannot locate the lsdaks executable after building.\n' >&2
    exit 1
fi
TABLE_DIR="${PROJECT_DIR}/data/tables/fortran_native"

MAX_ITER=5
FAILURES=0
CHECKS=0

pass() {
    CHECKS=$((CHECKS + 1))
    printf '  ok    %s\n' "$1"
}

fail() {
    CHECKS=$((CHECKS + 1))
    FAILURES=$((FAILURES + 1))
    printf '  FAIL  %s\n' "$1"
    if [ -n "${2:-}" ] && [ -f "$2" ]; then
        sed -n '1,200p' "$2" | sed 's/^/        | /'
    fi
}

# Rewrite an input file with a bounded iteration budget and an output prefix
# inside the work directory. Keys that are missing are appended (each of the
# four namelist groups is optional), keys that are present are overwritten -
# a second &scf group would be ignored, since the namelist reader takes the
# first match.
prepare_input() {
    local src="$1" dst="$2" prefix="$3"
    awk -v mi="${MAX_ITER}" -v pref="${prefix}" '
        /^[ \t]*&/ {
            grp = tolower($0); sub(/^[ \t]*&/, "", grp); sub(/[ \t].*$/, "", grp)
            if (grp == "scf")    seen_scf = 1
            if (grp == "output") seen_out = 1
            print; next
        }
        /^[ \t]*\/[ \t]*$/ {
            if (grp == "scf"    && !has_iter) print "  max_iter = " mi
            if (grp == "output" && !has_pref) print "  output_prefix = \047" pref "\047"
            grp = ""; print; next
        }
        tolower($0) ~ /^[ \t]*max_iter[ \t]*=/ {
            print "  max_iter = " mi; has_iter = 1; next
        }
        tolower($0) ~ /^[ \t]*output_prefix[ \t]*=/ {
            print "  output_prefix = \047" pref "\047"; has_pref = 1; next
        }
        { print }
        END {
            if (!seen_scf) print "&scf\n  max_iter = " mi "\n/"
            if (!seen_out) print "&output\n  output_prefix = \047" pref "\047\n/"
        }
    ' "${src}" > "${dst}"
}

run_lsdaks() {
    local input="$1" log="$2"
    (
        cd "${WORK_DIR}" || exit 1
        LSDAKS_TABLE_DIR="${TABLE_DIR}" "${EXECUTABLE}" --input "${input}"
    ) > "${log}" 2>&1
    return $?
}

shopt -s nullglob
EXAMPLES=("${PROJECT_DIR}"/examples/*.txt)
shopt -u nullglob

if [ ${#EXAMPLES[@]} -eq 0 ]; then
    fail "no example found in examples/"
fi

printf '\n== examples/*.txt must run AS COMMITTED (no preprocessing) ==\n'

# Pass 0: the bytes in the repository, untouched. No max_iter and no
# output_prefix override, because rewriting the file is exactly what used to
# hide a byte-level parser bug. The executable runs from WORK_DIR, so the
# committed relative prefixes stay isolated there.

for example in "${EXAMPLES[@]}"; do
    name="$(basename "${example}" .txt)"
    log="${WORK_DIR}/${name}.asis.log"

    run_lsdaks "${example}" "${log}"

    # As in the prepared pass below, the exit status is not the criterion: a
    # run that legitimately fails to converge exits 1.
    if ! grep -q "Starting Kohn-Sham SCF Cycle" "${log}"; then
        fail "${name} (as committed): did not reach the SCF cycle" "${log}"
    elif grep -q "^ *ERROR" "${log}"; then
        fail "${name} (as committed): reported an error before the SCF cycle" "${log}"
    else
        pass "${name} (as committed)"
    fi
done

printf '\n== examples/*.txt must run with a bounded iteration budget ==\n'

for example in "${EXAMPLES[@]}"; do
    name="$(basename "${example}" .txt)"
    prepared="${WORK_DIR}/${name}.txt"
    log="${WORK_DIR}/${name}.log"

    prepare_input "${example}" "${prepared}" "${WORK_DIR}/${name}"
    run_lsdaks "${prepared}" "${log}"

    # The exit status is NOT the criterion here: max_iter = 5 makes most of
    # these runs stop before convergence, which legitimately exits 1. What is
    # being checked is that the input was read, validated and turned into a
    # potential and a Hamiltonian, i.e. that the run reached the SCF cycle.
    if ! grep -q "Starting Kohn-Sham SCF Cycle" "${log}"; then
        fail "${name}: did not reach the SCF cycle" "${log}"
    elif grep -q "^ *ERROR" "${log}"; then
        fail "${name}: reported an error before the SCF cycle" "${log}"
    else
        pass "${name}"
    fi
done

printf '\n== a malformed input must be rejected ==\n'

# 2. Unknown key.
cat > "${WORK_DIR}/bad_key.txt" <<'EOF'
&system
  L = 10
  Nup = 5
  Ndown = 5
  U = 4.0
  no_such_key = 3
/
EOF
run_lsdaks "${WORK_DIR}/bad_key.txt" "${WORK_DIR}/bad_key.log"
status=$?
if [ "${status}" -eq 0 ]; then
    fail "an unknown namelist key must not be accepted" "${WORK_DIR}/bad_key.log"
elif ! grep -q "ERROR reading &system" "${WORK_DIR}/bad_key.log"; then
    fail "an unknown key must be reported by group and message" "${WORK_DIR}/bad_key.log"
elif ! grep -q "no_such_key" "${WORK_DIR}/bad_key.log"; then
    fail "the message must name the offending key" "${WORK_DIR}/bad_key.log"
else
    pass "unknown key rejected with a message naming it"
fi

# 3. Removed key: the message must say what to do, not just "cannot match".
cat > "${WORK_DIR}/obsolete_key.txt" <<'EOF'
&system
  L = 10
  Nup = 5
  Ndown = 5
  U = 4.0
/
&potential
  potential_type = 'random_uniform'
  disorder_strength = 1.0
  distribution = 'gaussian'
/
EOF
run_lsdaks "${WORK_DIR}/obsolete_key.txt" "${WORK_DIR}/obsolete_key.log"
status=$?
if [ "${status}" -eq 0 ]; then
    fail "the removed key 'distribution' must not be accepted" "${WORK_DIR}/obsolete_key.log"
elif ! grep -q "REMOVED" "${WORK_DIR}/obsolete_key.log"; then
    fail "'distribution' must be reported as removed, with the replacement" \
         "${WORK_DIR}/obsolete_key.log"
else
    pass "removed key 'distribution' rejected with a migration hint"
fi

# 4. Malformed value. gfortran reports this as plain end-of-file, exactly like
#    an absent group, so the parser has to tell the two apart by itself.
cat > "${WORK_DIR}/bad_value.txt" <<'EOF'
&system
  L = 10
  Nup = 5
  Ndown = 5
  U = 1.2.3
/
EOF
run_lsdaks "${WORK_DIR}/bad_value.txt" "${WORK_DIR}/bad_value.log"
status=$?
if [ "${status}" -eq 0 ]; then
    fail "a malformed value must not be accepted" "${WORK_DIR}/bad_value.log"
elif ! grep -q "ERROR reading &system" "${WORK_DIR}/bad_value.log"; then
    fail "a malformed value must be reported as a read error" "${WORK_DIR}/bad_value.log"
else
    pass "malformed value rejected with a message"
fi

# 5. Absent group: legitimate, and must stay legitimate.
cat > "${WORK_DIR}/no_groups.txt" <<EOF
&system
  L = 10
  Nup = 5
  Ndown = 5
  U = 4.0
/
&scf
  max_iter = ${MAX_ITER}
  verbose = .false.
/
&output
  output_prefix = '${WORK_DIR}/no_groups'
/
EOF
run_lsdaks "${WORK_DIR}/no_groups.txt" "${WORK_DIR}/no_groups.log"
if ! grep -q "Starting Kohn-Sham SCF Cycle" "${WORK_DIR}/no_groups.log"; then
    fail "an absent &potential group must be accepted" "${WORK_DIR}/no_groups.log"
elif ! grep -q "group &potential not found" "${WORK_DIR}/no_groups.log"; then
    fail "an absent group must be reported as a note" "${WORK_DIR}/no_groups.log"
else
    pass "absent &potential group accepted, with a note"
fi

printf '\n== a converged run whose output is lost must fail ==\n'

# 6. The SCF converges and every output file fails to be written. The numbers
#    are gone, so only the exit status can tell the caller.
cat > "${WORK_DIR}/unwritable.txt" <<'EOF'
&system
  L = 10
  Nup = 5
  Ndown = 5
  U = 4.0
/
&scf
  max_iter = 200
  verbose = .false.
/
&output
  output_prefix = '/lsdaks_no_such_directory/out'
/
EOF
run_lsdaks "${WORK_DIR}/unwritable.txt" "${WORK_DIR}/unwritable.log"
status=$?
# The guard below must NOT be `grep -q "CONVERGED"`: that pattern is a
# substring of "NOT CONVERGED", so it is satisfied by the very run it is meant
# to exclude (output_writer.f90 warns about exactly this). The status line is
# matched with its ✓ marker instead, which only the converged branch prints.
if ! grep -q "Status:.*✓ CONVERGED" "${WORK_DIR}/unwritable.log"; then
    fail "the write-failure check needs a run that converges" "${WORK_DIR}/unwritable.log"
elif [ "${status}" -eq 0 ]; then
    fail "a converged run whose output could not be written must exit non-zero" \
         "${WORK_DIR}/unwritable.log"
elif ! grep -q "could NOT be written" "${WORK_DIR}/unwritable.log"; then
    fail "the lost output must be reported explicitly" "${WORK_DIR}/unwritable.log"
else
    pass "converged run with unwritable output exits non-zero"
fi

printf '\n== the table generator must not persist an invalid table ==\n'

# 7. Invariant: (exit 0 and a finite table) or (exit non-zero and no file).
#    The generator's grid is fixed (50 x 51, U=4 takes a few seconds), so the
#    check runs it for real instead of relying on a unit-level stub.
GEN_REL="$(fpm run generate_xc_table --runner echo | tail -n 1)"
case "${GEN_REL}" in
    /*) GENERATOR="${GEN_REL}" ;;
    *)  GENERATOR="${PROJECT_DIR}/${GEN_REL}" ;;
esac
GEN_DIR="${WORK_DIR}/gen"
mkdir -p "${GEN_DIR}"
if [ ! -x "${GENERATOR}" ]; then
    fail "cannot locate the generate_xc_table executable"
else
    "${GENERATOR}" --U 4.0 --output "${GEN_DIR}" > "${WORK_DIR}/generator.log" 2>&1
    status=$?
    GEN_FILE="${GEN_DIR}/xc_table_u4.00.dat"
    if [ "${status}" -eq 0 ]; then
        if [ ! -f "${GEN_FILE}" ]; then
            fail "generator exited 0 but wrote no table" "${WORK_DIR}/generator.log"
        # The table is an unformatted stream (see table_io::write_fortran_table):
        # U, the two grid sizes, n_grid, m_grid, then exc, vxc_up, vxc_down as
        # raw little-endian doubles. Decode the three data arrays and look for
        # NaN or ±Inf.
        elif python3 - "${GEN_FILE}" <<'PY'
import struct, sys
data = open(sys.argv[1], "rb").read()
# header: U (f8), n_points_n, n_points_m (i4, i4)
nn, nm = struct.unpack_from("<ii", data, 8)
off = 16 + 8 * nn + 8 * nn * nm            # skip n_grid and m_grid
vals = struct.unpack_from("<%dd" % (3 * nn * nm), data, off)
bad = sum(1 for v in vals if v != v or v in (float("inf"), float("-inf")))
sys.exit(1 if bad else 0)
PY
        then
            pass "generator exited 0 and the table it wrote is finite"
        else
            fail "generator exited 0 but the table contains NaN/Inf" "${WORK_DIR}/generator.log"
        fi
    else
        if [ -e "${GEN_FILE}" ]; then
            fail "generator exited non-zero but left a table behind" "${WORK_DIR}/generator.log"
        elif ! grep -q "non-finite" "${WORK_DIR}/generator.log"; then
            fail "generator failure must be reported as non-finite entries" "${WORK_DIR}/generator.log"
        else
            pass "generator exited non-zero with a non-finite table and wrote nothing"
        fi
    fi
fi

printf '\n%d checks, %d failures\n\n' "${CHECKS}" "${FAILURES}"

if [ "${FAILURES}" -ne 0 ]; then
    exit 1
fi
exit 0
