#!/bin/bash
#SBATCH -A naiss2025-5-732
#SBATCH -t 00-06:00:00
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH -J rev7_rare_common
#SBATCH -o rev7_rare_common.%j.out
#SBATCH -e rev7_rare_common.%j.err

python evaluate_rare_common.py
