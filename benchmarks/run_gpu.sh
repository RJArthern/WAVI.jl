#!/bin/bash -l
# WAVI single-GPU GPUSpec run (not a scaling sweep).
#
# CPU thread and MPI sweeps are in benchmarks/run_scaling.sh.
# This script is configured for the BAS internal HPC:
# one NVIDIA GPU on bsl-node-s22 (8 x V100S). GPUSpec currently only uses one device;
# do not set CUDA_VISIBLE_DEVICES when Slurm has already allocated `--gres=gpu:1`.
#
# Usage (login node): sbatch benchmarks/run_gpu.sh
# Usage (already on the GPU node): bash benchmarks/run_gpu.sh
# Optional args: tag prefix, then driver name.
#
# Default node is bsl-node-s22 (8 x V100S). Override on the sbatch line:
#   sbatch --nodelist=bsl-node-s20 benchmarks/run_gpu.sh   # Currently has 2 x A40 GPUs
#   sbatch --nodelist=bsl-node-s21 benchmarks/run_gpu.sh   # Currently has 2 x T4 GPUs

#SBATCH -A gpu
#SBATCH -p gpu
#SBATCH -J wavi_gpu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --nodelist=bsl-node-s22
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00
#SBATCH --output=benchmarks/output/ismip7_16km_synthetic/slurm-gpu-%j.out

# GPUSpec itself is one process on one device (no mpiexec yet).
# Do not `module load cuda`: CUDA.jl vendors its own toolkit and so doing this would
# cause CUDA.jl to warn about mixing ABIs (since it adds the site module's libraries
# to LD_LIBRARY_PATH). Just the NVIDIA driver on the GPU node is enough.

set -e

TAG_PREFIX=${1:-"ka.jl_gpuspec_v1"}
DRIVER=${2:-"ismip7_16km_synthetic"}
NODE_TAG=$(hostname -s | sed 's/^bsl-node-//')
RUN_TAG="${TAG_PREFIX}_${NODE_TAG}_gpu1"

mkdir -p "benchmarks/output/${DRIVER}"

echo "=========================================="
echo "Starting WAVI GPUSpec run"
echo "Node: $(hostname)"
echo "Commit: $(git rev-parse HEAD)"
echo "Driver: $DRIVER"
echo "Tag Prefix: $TAG_PREFIX"
echo "Run tag: $RUN_TAG"
echo "CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES:-<unset>}"
echo "=========================================="

if command -v nvidia-smi >/dev/null 2>&1; then
    nvidia-smi
else
    echo "nvidia-smi not found on PATH"
fi

# Making sure CUDA.jl is in the benchmarks project (Note: not in WAVI [deps]).
julia --project=benchmarks -e '
using Pkg
if !haskey(Pkg.project().dependencies, "CUDA")
    println("Adding CUDA.jl to the benchmarks project...")
    Pkg.add("CUDA")
end
'

echo "Running GPUSpec (1 GPU)..."
julia --project=benchmarks -t 1 benchmarks/run.jl run gpu "$DRIVER" \
    --tag "$RUN_TAG" --output-group "$TAG_PREFIX"

echo "GPUSpec run complete!"
