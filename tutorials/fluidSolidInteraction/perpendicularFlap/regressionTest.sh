#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ============================================================
# perpendicularFlap FSI regression test
# ============================================================

# ------------------------------------------------------------
# Regression tolerances
# ------------------------------------------------------------

DISP_MAX_TOL=1e-6      # max displacement absolute tolerance

# Reference values
REF_MAX_DISP=0.118329

# Log files
ALLRUN_LOGFILE="log.Allrun"

# Data files
DISP_FILE="postProcessing/0/solidPointDisplacement_displacement.dat"

echo "============================================================"
echo "perpendicularFlap FSI regression test"
echo "Max displacement difference < ${DISP_MAX_TOL}"
echo "============================================================"
echo

# ------------------------------------------------------------
# Clean & run case
# ------------------------------------------------------------

./Allclean > /dev/null 2>&1 || true
./Allrun > "${ALLRUN_LOGFILE}" 2>&1

# ------------------------------------------------------------
# Extract helpers
# ------------------------------------------------------------

extract_max_displacement() {
    awk '{print $2}' "${DISP_FILE}" | sort -g | tail -1
}

abs() {
    awk -v x="$1" 'BEGIN {print (x < 0 ? -x : x)}'
}

# ------------------------------------------------------------
# Extract values
# ------------------------------------------------------------

max_disp=$(extract_max_displacement)

if [[ -z "${max_disp}" ]]; then
    echo "FAIL: Could not extract regression quantities"
    exit 1
fi

disp_diff=$(awk "BEGIN {print ${max_disp} - ${REF_MAX_DISP}}")
disp_diff_abs=$(abs "${disp_diff}")

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

failures=0

if awk "BEGIN {exit !(${disp_diff_abs} < ${DISP_MAX_TOL})}"; then
    printf "PASS: max displacement = %.6g (Δ = %.3g)\n" \
        "${max_disp}" "${disp_diff_abs}"
else
    printf "FAIL: max displacement = %.6g (Δ = %.3g)\n" \
        "${max_disp}" "${disp_diff_abs}"
    failures=$((failures + 1))
fi

# Clean case again
./Allclean > /dev/null 2>&1 || true

echo
if (( failures == 0 )); then
    echo "============================================================"
    echo "Regression test PASSED"
    echo "============================================================"
    exit 0
else
    echo "============================================================"
    echo "Regression test FAILED (${failures} checks)"
    echo "============================================================"
    exit 1
fi
