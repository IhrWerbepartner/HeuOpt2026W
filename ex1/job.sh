#!/bin/bash
#SBATCH --job-name=rckst
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x-%j.out

module load julia
julia --project=. src/main.jl "$@"