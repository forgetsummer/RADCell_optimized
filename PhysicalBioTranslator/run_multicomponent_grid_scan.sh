#!/bin/bash
# Grid scan over (Nc, omega_p, lambda_p) for multi-component energy model.
# Runs test_multicomponent_validation for each combination and collects results.

set -e

# Build directory -- adjust if different
BUILD_DIR="${BUILD_DIR:-/home/user/Geant4Projects/CellStateTransitionPaper/RAD-build}"
PIDE_DIR="${PIDE_DIR:-/home/user/Geant4Projects/CellStateTransitionPaper/RADCellSimulation/PIDE3.4}"
BASE_OUTPUT="${BUILD_DIR}/PhysicalBioTranslator/multicomponent_grid_scan"
EXECUTABLE="${BUILD_DIR}/PhysicalBioTranslator/test_multicomponent_validation"

CELLS=200
REPLICATES=3
MAX_CELL_TYPES=10

# Grid values
NC_VALUES=(15 30 60)
OMEGA_P_VALUES=(1.0 1.5 2.0 3.0)
LAMBDA_P_VALUES=(0.02 0.04 0.07)

# Also run baseline (Nc=99999 effectively disables Er/Ep split since h(N)~0)
RUN_BASELINE=true

mkdir -p "$BASE_OUTPUT"

echo "============================================"
echo "  Multi-Component Energy Grid Scan"
echo "============================================"
echo "Build dir:      $BUILD_DIR"
echo "PIDE dir:       $PIDE_DIR"
echo "Output base:    $BASE_OUTPUT"
echo "Grid:           ${#NC_VALUES[@]} x ${#OMEGA_P_VALUES[@]} x ${#LAMBDA_P_VALUES[@]} = $(( ${#NC_VALUES[@]} * ${#OMEGA_P_VALUES[@]} * ${#LAMBDA_P_VALUES[@]} )) combinations"
echo "Cells/sim:      $CELLS"
echo "Replicates:     $REPLICATES"
echo ""

if [ ! -f "$EXECUTABLE" ]; then
    echo "ERROR: Executable not found: $EXECUTABLE"
    echo "Please build first: cd $BUILD_DIR && make test_multicomponent_validation"
    exit 1
fi

COMBINED_CSV="$BASE_OUTPUT/grid_scan_results.csv"
echo "Nc,omega_p,lambda_p,CellLine,Alpha_LQ,Beta_LQ,AvgLog10Error,MaxLog10Error" > "$COMBINED_CSV"

RUN_COUNT=0
TOTAL_RUNS=$(( ${#NC_VALUES[@]} * ${#OMEGA_P_VALUES[@]} * ${#LAMBDA_P_VALUES[@]} ))

if [ "$RUN_BASELINE" = true ]; then
    TOTAL_RUNS=$(( TOTAL_RUNS + 1 ))
    echo "--- Running baseline (Nc=99999, omega_p=1.0, lambda_p=0.0) ---"
    OUTDIR="$BASE_OUTPUT/baseline"
    mkdir -p "$OUTDIR"
    "$EXECUTABLE" \
        --pide-dir "$PIDE_DIR" \
        --output-dir "$OUTDIR" \
        --cells "$CELLS" \
        --replicates "$REPLICATES" \
        --max-cell-types "$MAX_CELL_TYPES" \
        --nc 99999 \
        --omega-p 1.0 \
        --lambda-p 0.001 \
        2>&1 | tee "$OUTDIR/log.txt"

    if [ -f "$OUTDIR/mc_validation_summary.csv" ]; then
        tail -n +2 "$OUTDIR/mc_validation_summary.csv" | while IFS=',' read -r CellLine PhotonRad Alpha Beta ABR e2 e3 a kerr T21 T23 Conv AvgErr MaxErr Nc Op Lp; do
            echo "99999,1.0,0.001,$CellLine,$Alpha,$Beta,$AvgErr,$MaxErr" >> "$COMBINED_CSV"
        done
    fi
    RUN_COUNT=$((RUN_COUNT + 1))
    echo ""
fi

for NC in "${NC_VALUES[@]}"; do
    for OP in "${OMEGA_P_VALUES[@]}"; do
        for LP in "${LAMBDA_P_VALUES[@]}"; do
            RUN_COUNT=$((RUN_COUNT + 1))
            echo "--- Run $RUN_COUNT/$TOTAL_RUNS: Nc=$NC omega_p=$OP lambda_p=$LP ---"

            OUTDIR="$BASE_OUTPUT/nc${NC}_op${OP}_lp${LP}"
            mkdir -p "$OUTDIR"

            "$EXECUTABLE" \
                --pide-dir "$PIDE_DIR" \
                --output-dir "$OUTDIR" \
                --cells "$CELLS" \
                --replicates "$REPLICATES" \
                --max-cell-types "$MAX_CELL_TYPES" \
                --nc "$NC" \
                --omega-p "$OP" \
                --lambda-p "$LP" \
                2>&1 | tee "$OUTDIR/log.txt"

            if [ -f "$OUTDIR/mc_validation_summary.csv" ]; then
                tail -n +2 "$OUTDIR/mc_validation_summary.csv" | while IFS=',' read -r CellLine PhotonRad Alpha Beta ABR e2 e3 a kerr T21 T23 Conv AvgErr MaxErr Nc_csv Op_csv Lp_csv; do
                    echo "$NC,$OP,$LP,$CellLine,$Alpha,$Beta,$AvgErr,$MaxErr" >> "$COMBINED_CSV"
                done
            fi

            echo ""
        done
    done
done

echo "============================================"
echo "  Grid Scan Complete"
echo "============================================"
echo "Total runs: $RUN_COUNT"
echo "Combined results: $COMBINED_CSV"
echo ""
echo "Next step: python3 analyze_multicomponent_results.py $BASE_OUTPUT"
