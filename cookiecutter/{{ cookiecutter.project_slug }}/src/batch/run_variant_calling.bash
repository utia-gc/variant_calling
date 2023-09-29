#!/bin/bash
#SBATCH --job-name={{ cookiecutter.project_slug }}_variant-calling
#SBATCH --error=job_logs/%x.e%j
#SBATCH --output=job_logs/%x.o%j
#SBATCH --mail-user={{ cookiecutter.email }}
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --account={{ cookiecutter.slurm_account }}
#SBATCH --partition={{ cookiecutter.slurm_partition }}
#SBATCH --qos={{ cookiecutter.slurm_qos }}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01-00:00:00


# set nextflow options for execution via slurm
export NXF_OPTS="-Xms500M -Xmx2G"
export NXF_ANSI_LOG=false

# run pipeline
nextflow run utia-gc/variant_calling \
    -latest \
    -revision v0.0.0.9001 \
    -profile isaac_tff,exploratory
