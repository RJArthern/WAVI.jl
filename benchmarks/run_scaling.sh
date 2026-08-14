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
#SBATCH --time=03:00:00

# Ensure we use the correct MPI exec wrapper for the BAS HPC environment
MPIEXEC=mpiexecjl

echo "=========================================="
echo "Starting WAVI Scaling Benchmark Suite"
echo "Node: $(hostname)"
echo "Commit: $(git rev-parse HEAD)"
echo "=========================================="

echo "Running 1-core MPI baseline..."
$MPIEXEC -n 1 julia --project=benchmarks -t 1 benchmarks/run.jl run mpi ismip7_16km_synthetic --px 1 --py 1 --tag "baseline_sparse_1core"

echo "Running 2-core (2x1) scaling..."
$MPIEXEC -n 2 julia --project=benchmarks -t 1 benchmarks/run.jl run mpi ismip7_16km_synthetic --px 2 --py 1 --tag "baseline_sparse_2core"

echo "Running 4-core (2x2) scaling..."
$MPIEXEC -n 4 julia --project=benchmarks -t 1 benchmarks/run.jl run mpi ismip7_16km_synthetic --px 2 --py 2 --tag "baseline_sparse_4core"

echo "Running 8-core (4x2) scaling..."
$MPIEXEC -n 8 julia --project=benchmarks -t 1 benchmarks/run.jl run mpi ismip7_16km_synthetic --px 4 --py 2 --tag "baseline_sparse_8core"

echo "Running 16-core (4x4) scaling..."
$MPIEXEC -n 16 julia --project=benchmarks -t 1 benchmarks/run.jl run mpi ismip7_16km_synthetic --px 4 --py 4 --tag "baseline_sparse_16core"

echo "Running 36-core (6x6) scaling..."
$MPIEXEC -n 36 julia --project=benchmarks -t 1 benchmarks/run.jl run mpi ismip7_16km_synthetic --px 6 --py 6 --tag "baseline_sparse_36core"

echo "Scaling benchmarks complete!"
