#!/bin/bash
set -e

# Script location
S_DIR="$(cd "$(dirname "$0")" && pwd)"

# Base directory for all test runs
BASE_DIR="$S_DIR/outputs"
mkdir -p "$BASE_DIR"

echo "--------------------------------------------------"
echo "Running MISMIP_PLUS - BasicSpec (Serial)"
echo "--------------------------------------------------"
mkdir -p "$BASE_DIR/basic"
cd "$BASE_DIR/basic"
julia --project="$S_DIR" "$S_DIR/MISMIP_PLUS.jl"
mv outputs/* .
rmdir outputs

echo "--------------------------------------------------"
echo "Running MISMIP_PLUS - ThreadedSpec"
echo "--------------------------------------------------"
for t in 2 3 5 7; do
    echo "  Threads: $t"
    mkdir -p "$BASE_DIR/threaded_t$t"
    cd "$BASE_DIR/threaded_t$t"
    julia --project="$S_DIR" -t "$t" "$S_DIR/MISMIP_PLUS.jl"
    mv outputs/* .
    rmdir outputs
done

echo "--------------------------------------------------"
echo "Running MISMIP_PLUS - MPISpec"
echo "--------------------------------------------------"
for n in 2 3 5 7; do
    echo "  MPI Processes: $n"
    mkdir -p "$BASE_DIR/mpi_n$n"
    cd "$BASE_DIR/mpi_n$n"
    mpiexecjl -n "$n" julia --project="$S_DIR" "$S_DIR/MISMIP_PLUS.jl"
    mv outputs/* .
    rmdir outputs
done

echo "--------------------------------------------------"
echo "Running MISMIP_PLUS - GPUSpec"
echo "--------------------------------------------------"
mkdir -p "$BASE_DIR/gpu"
cd "$BASE_DIR/gpu"
WAVI_USE_GPU=1 julia --project="$S_DIR" -t 1 "$S_DIR/MISMIP_PLUS.jl"
if [ -d outputs ]; then
    mv outputs/* .
    rmdir outputs
fi

echo "--------------------------------------------------"
echo "All MISMIP+ tests completed. Outputs in $BASE_DIR"
echo "--------------------------------------------------"
