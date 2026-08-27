#!/bin/bash -l
# WAVI MPISpec + GPU sweep: one GPU per MPI rank (not a CPU thread sweep).
#
# CPU thread and CPU MPI sweeps are in benchmarks/run_scaling.sh.
# Single-GPU GPUSpec is in benchmarks/run_gpu.sh.
# Default node is bsl-node-s22 (8 x V100S): request 4 GPUs and sweep 1-4 ranks.
# bsl-node-s20 (2 x A40) and bsl-node-s21 (2 x T4) only have two GPUs: override
# --gres and --ntasks. The sweep is 1..N for N GPUs allocated by Slurm.
#
# Usage (login node): sbatch benchmarks/run_mpi_gpu.sh
# Usage (already on a GPU node): bash benchmarks/run_mpi_gpu.sh
# Optional args: tag prefix, then driver name.
# GPU count is SLURM_JOB_GRES when set, otherwise nvidia-smi.
#
#   sbatch --nodelist=bsl-node-s20 --gres=gpu:2 --ntasks=2 benchmarks/run_mpi_gpu.sh
#   sbatch --nodelist=bsl-node-s21 --gres=gpu:2 --ntasks=2 benchmarks/run_mpi_gpu.sh

#SBATCH -A gpu
#SBATCH -p gpu
#SBATCH -J wavi_mpi_gpu
#SBATCH -N 1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --nodelist=bsl-node-s22
#SBATCH --gres=gpu:4
#SBATCH --time=04:00:00
#SBATCH --output=benchmarks/output/ismip7_16km_synthetic/slurm-mpi-gpu-%j.out

set -e

TAG_PREFIX=${1:-"ka.jl_haar_rap_pack_v7"}
DRIVER=${2:-"ismip7_16km_synthetic"}
NODE_TAG=$(hostname -s | sed 's/^bsl-node-//')
MPIEXEC=mpiexecjl

# Allocated GPUs: Slurm --gres if present, else nvidia-smi (plain bash on a GPU system).
if [[ -n "${SLURM_JOB_GRES:-}" ]]; then
    MAX_GPUS=${SLURM_JOB_GRES##*:}
    MAX_GPUS=${MAX_GPUS%%(*}
elif command -v nvidia-smi >/dev/null 2>&1; then
    MAX_GPUS=$(nvidia-smi -L 2>/dev/null | grep -c '^GPU' || true)
else
    MAX_GPUS=0
fi
if [[ -z "$MAX_GPUS" || "$MAX_GPUS" -lt 1 ]]; then
    echo "No GPU count (SLURM_JOB_GRES='${SLURM_JOB_GRES:-}', nvidia-smi missing or empty). Submit with sbatch --gres=gpu:N, or run on a node with nvidia-smi." >&2
    exit 1
fi

mkdir -p "benchmarks/output/${DRIVER}"

echo "=========================================="
echo "Starting WAVI MPISpec GPU sweep"
echo "Node: $(hostname)"
echo "Commit: $(git rev-parse HEAD)"
echo "Driver: $DRIVER"
echo "Tag Prefix: $TAG_PREFIX"
echo "Sweep: 1..${MAX_GPUS} ranks (one GPU per rank)"
echo "SLURM_JOB_GRES: ${SLURM_JOB_GRES:-<unset>}"
echo "CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES:-<unset>}"
echo "=========================================="

if command -v nvidia-smi >/dev/null 2>&1; then
    nvidia-smi
else
    echo "nvidia-smi not found on PATH"
fi

julia --project=benchmarks -e '
using Pkg
if !haskey(Pkg.project().dependencies, "CUDA")
    println("Adding CUDA.jl to the benchmarks project...")
    Pkg.add("CUDA")
end
'

run_n() {
    local n="$1"
    local run_tag="${TAG_PREFIX}_${NODE_TAG}_gpu${n}"
    echo "Running MPISpec GPU (${n} ranks, ${n}x1)..."
    $MPIEXEC -n "$n" julia --project=benchmarks -t 1 benchmarks/run.jl run mpi_gpu "$DRIVER" \
        --px "$n" --py 1 --tag "$run_tag" --output-group "$TAG_PREFIX"
}

for n in $(seq 1 "$MAX_GPUS"); do
    run_n "$n"
done

echo "MPISpec GPU sweep complete!"
