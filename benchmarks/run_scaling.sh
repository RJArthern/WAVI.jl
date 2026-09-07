#!/bin/bash
# WAVI Scaling Benchmark Sweep
#
# This script is configured for the British Antarctic Survey (BAS) HPC.
# It assumes submission to the 'medium' partition on a node with 36
# physical cores.
#
# Usage (Interactive node): bash benchmarks/run_scaling.sh
# Usage (Login node): sbatch benchmarks/run_scaling.sh

#SBATCH -A medium
#SBATCH -p medium
#SBATCH -J wavi_scale
#SBATCH -N 1
#SBATCH -n 36
#SBATCH --time=08:00:00
#SBATCH --output=benchmarks/output/ismip7_16km_synthetic/slurm-%j.out

# Ensure we use the correct MPI exec wrapper for the BAS HPC environment
MPIEXEC=mpiexecjl

# Accept tag prefix from command line argument 1
TAG_PREFIX=${1:-"ka.jl_mpi_schwarz_v15"}

# Accept driver from command line argument 2
DRIVER=${2:-"ismip7_16km_synthetic"}

echo "=========================================="
echo "Starting WAVI Scaling Benchmark Suite"
echo "Node: $(hostname)"
echo "Commit: $(git rev-parse HEAD)"
echo "Driver: $DRIVER"
echo "Tag Prefix: $TAG_PREFIX"
echo "=========================================="

echo "Running BasicSpec (serial) baseline..."
julia --project=benchmarks -t 1 benchmarks/run.jl run basic $DRIVER --tag "${TAG_PREFIX}_serial" --output-group "$TAG_PREFIX"

echo "Running BasicSpec 2-thread scaling..."
julia --project=benchmarks -t 2 benchmarks/run.jl run basic $DRIVER --tag "${TAG_PREFIX}_threads_2" --output-group "$TAG_PREFIX"

echo "Running BasicSpec 4-thread scaling..."
julia --project=benchmarks -t 4 benchmarks/run.jl run basic $DRIVER --tag "${TAG_PREFIX}_threads_4" --output-group "$TAG_PREFIX"

echo "Running BasicSpec 8-thread scaling..."
julia --project=benchmarks -t 8 benchmarks/run.jl run basic $DRIVER --tag "${TAG_PREFIX}_threads_8" --output-group "$TAG_PREFIX"

echo "Running BasicSpec 16-thread scaling..."
julia --project=benchmarks -t 16 benchmarks/run.jl run basic $DRIVER --tag "${TAG_PREFIX}_threads_16" --output-group "$TAG_PREFIX"

echo "Running BasicSpec 36-thread scaling..."
julia --project=benchmarks -t 36 benchmarks/run.jl run basic $DRIVER --tag "${TAG_PREFIX}_threads_36" --output-group "$TAG_PREFIX"

echo "Running 1-core MPI baseline..."
$MPIEXEC -n 1 julia --project=benchmarks -t 1 benchmarks/run.jl run mpi $DRIVER --px 1 --py 1 --tag "${TAG_PREFIX}_1core" --output-group "$TAG_PREFIX"

echo "Running 2-core (2x1) scaling..."
$MPIEXEC -n 2 julia --project=benchmarks -t 1 benchmarks/run.jl run mpi $DRIVER --px 2 --py 1 --tag "${TAG_PREFIX}_2core" --output-group "$TAG_PREFIX"

echo "Running 4-core (2x2) scaling..."
$MPIEXEC -n 4 julia --project=benchmarks -t 1 benchmarks/run.jl run mpi $DRIVER --px 2 --py 2 --tag "${TAG_PREFIX}_4core" --output-group "$TAG_PREFIX"

echo "Running 8-core (4x2) scaling..."
$MPIEXEC -n 8 julia --project=benchmarks -t 1 benchmarks/run.jl run mpi $DRIVER --px 4 --py 2 --tag "${TAG_PREFIX}_8core" --output-group "$TAG_PREFIX"

echo "Running 16-core (4x4) scaling..."
$MPIEXEC -n 16 julia --project=benchmarks -t 1 benchmarks/run.jl run mpi $DRIVER --px 4 --py 4 --tag "${TAG_PREFIX}_16core" --output-group "$TAG_PREFIX"

echo "Running 36-core (6x6) scaling..."
$MPIEXEC -n 36 julia --project=benchmarks -t 1 benchmarks/run.jl run mpi $DRIVER --px 6 --py 6 --tag "${TAG_PREFIX}_36core" --output-group "$TAG_PREFIX"

echo "Scaling benchmarks complete!"
