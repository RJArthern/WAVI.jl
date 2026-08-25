#!/bin/bash
set -e

# Script location
S_DIR="$(cd "$(dirname "$0")" && pwd)"

# Base directory for all test runs
BASE_DIR="$S_DIR/outputs"
mkdir -p "$BASE_DIR"

# Move Julia's outputs/<kind>/ tree into the current run directory.
# Safe to re-run: a leftover gpu/ or mpi/ folder is replaced.
flatten_run_outputs() {
    [[ -d outputs ]] || return 0
    local item dest
    for item in outputs/*; do
        [[ -e "$item" ]] || return 0
        dest="./$(basename "$item")"
        [[ -e "$dest" ]] && rm -rf "$dest"
        mv "$item" "$dest"
    done
    rmdir outputs
}

echo "--------------------------------------------------"
echo "Running MISMIP_PLUS - BasicSpec (Serial)"
echo "--------------------------------------------------"
mkdir -p "$BASE_DIR/basic"
cd "$BASE_DIR/basic"
julia --project="$S_DIR" "$S_DIR/MISMIP_PLUS.jl"
flatten_run_outputs

echo "--------------------------------------------------"
echo "Running MISMIP_PLUS - ThreadedSpec"
echo "--------------------------------------------------"
for t in 2 3 5 7; do
    echo "  Threads: $t"
    mkdir -p "$BASE_DIR/threaded_t$t"
    cd "$BASE_DIR/threaded_t$t"
    julia --project="$S_DIR" -t "$t" "$S_DIR/MISMIP_PLUS.jl"
    flatten_run_outputs
done

echo "--------------------------------------------------"
echo "Running MISMIP_PLUS - MPISpec"
echo "--------------------------------------------------"
for n in 2 3 5 7; do
    echo "  MPI Processes: $n"
    mkdir -p "$BASE_DIR/mpi_n$n"
    cd "$BASE_DIR/mpi_n$n"
    mpiexecjl -n "$n" julia --project="$S_DIR" "$S_DIR/MISMIP_PLUS.jl"
    flatten_run_outputs
done

echo "--------------------------------------------------"
echo "Running MISMIP_PLUS - GPUSpec"
echo "--------------------------------------------------"
mkdir -p "$BASE_DIR/gpu"
cd "$BASE_DIR/gpu"
WAVI_USE_GPU=1 julia --project="$S_DIR" -t 1 "$S_DIR/MISMIP_PLUS.jl"
flatten_run_outputs

# MPISpec + one GPU per rank, only if this allocation has multiple GPUs.
# Run as one Slurm task; this script starts mpiexecjl itself.
# `cpus-per-task` covering the largest CPUs needed (ThreadedSpec -t 7, MPISpec -n 7):
# srun -A gpu -p gpu --ntasks=1 --cpus-per-task=8 --nodelist=bsl-node-s22 \
#      --mem=128G --gres=gpu:4 --time=04:00:00 ./run_mismip_plus.sh
if [[ -n "${SLURM_JOB_GRES:-}" ]]; then
    ngpus=${SLURM_JOB_GRES##*:}
    ngpus=${ngpus%%(*}
elif command -v nvidia-smi >/dev/null 2>&1; then
    ngpus=$(nvidia-smi -L 2>/dev/null | grep -c '^GPU' || true)
else
    ngpus=0
fi

if [[ "${ngpus:-0}" -gt 1 ]]; then
    echo "--------------------------------------------------"
    echo "Running MISMIP_PLUS - MPISpec + GPU ($ngpus devices)"
    echo "--------------------------------------------------"
    echo "  Skipping 1 rank (that is GPUSpec, already run above). MPISpec+GPU starts at 2 ranks."
    for n in $(seq 2 "$ngpus"); do
        echo "  MPISpec + GPU: $n ranks / GPUs"
        mkdir -p "$BASE_DIR/mpi_gpu_n$n"
        cd "$BASE_DIR/mpi_gpu_n$n"
        WAVI_USE_GPU=1 mpiexecjl -n "$n" julia --project="$S_DIR" -t 1 "$S_DIR/MISMIP_PLUS.jl"
        flatten_run_outputs
    done
else
    echo "Skipping MPI+GPU sweep (need more than one GPU; found ${ngpus:-0})"
fi

echo "--------------------------------------------------"
echo "All MISMIP+ tests completed. Outputs in $BASE_DIR"
echo "--------------------------------------------------"
