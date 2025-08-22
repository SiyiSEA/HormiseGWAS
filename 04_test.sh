#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH -A Research_Project-MRC190311 # research project to submit under.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --output=Jobreports/mungeT2D.out
#SBATCH --error=Jobreports/mungeT2D.err
#SBATCH --job-name=mungeT2D
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --ntasks=2
#SBATCH --time=0-15:00:00

RscriptsPath="/lustre/home/sww208/GoDMC/HormiseGWAS/Rscripts"
DataPath="/lustre/home/sww208/GoDMC/GADatasets/GWASstatsFiltered"

module purge
source "/gpfs/ts0/shared/software/Miniconda3/23.5.2-0/etc/profile.d/conda.sh"
conda activate RGreatSquareRoot

Rscript $RscriptsPath/munge_MSS.R $DataPath T2D_2024.b37.filtered.SNP T2D_2024