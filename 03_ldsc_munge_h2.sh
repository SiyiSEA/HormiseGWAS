#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH -A Research_Project-MRC190311 # research project to submit under.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --output=Jobreports/ldscMugne.out
#SBATCH --error=Jobreports/ldscMugne.err
#SBATCH --job-name=ldscMugne
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --ntasks=2
#SBATCH --time=0-10:00:00
#SBATCH --array=2



#define the path
w_hm3="/lustre/home/sww208/GoDMC/GADatasets/CorrectGWAS10/Resources/w_hm3.snplist"
ref_ld="/lustre/home/sww208/GoDMC/GADatasets/CorrectGWAS10/Resources/baselineLD_v2.2/baselineLD."
w_ld="/lustre/home/sww208/GoDMC/GADatasets/CorrectGWAS10/Resources/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."

source "/gpfs/ts0/shared/software/Miniconda3/23.5.2-0/etc/profile.d/conda.sh"
conda activate ldsc
LDSC="/lustre/home/sww208/Software/ldsc"


if [ "${SLURM_ARRAY_TASK_ID}" == 1 ] ;then
    DataPath="/lustre/home/sww208/GoDMC/GADatasets/OtherGWASstats"
    echo "Dealing with datasets in $DataPath"
elif [ "${SLURM_ARRAY_TASK_ID}" == 2 ] ;then
    DataPath="/lustre/home/sww208/GoDMC/GADatasets/GWASstatsFiltered"
    echo "Dealing with datasets in $DataPath"
fi

cd ${DataPath} || exit

ldsc_munge() {
    local ID=$1
    local HasNcol=$2
    local N=$3

    echo "Munging ${ID}..."

    if [ $HasNcol = true ]; then
      ${LDSC}/munge_sumstats.py \
          --sumstats ${ID}_MSS.tsv \
          --snp SNP \
          --N-col $N \
          --a1 A1 \
          --a2 A2 \
          --frq FRQ \
          --p P \
          --chunksize 50000 \
          --signed-sumstats Z,0 \
          --merge-alleles ${w_hm3} \
          --out ${ID}
          
    elif [ $HasNcol = false ]; then
      ${LDSC}/munge_sumstats.py \
          --sumstats ${ID}_MSS.tsv \
          --snp SNP \
          --N $N \
          --a1 A1 \
          --a2 A2 \
          --frq FRQ \
          --p P \
          --chunksize 50000 \
          --signed-sumstats Z,0 \
          --merge-alleles ${w_hm3} \
          --out ${ID}

    elif [ $HasNcol = case ]; then
      ${LDSC}/munge_sumstats.py \
          --sumstats ${ID}_MSS.tsv \
          --snp SNP \
          --N-cas-col $N \
          --N-con-col $4 \
          --a1 A1 \
          --a2 A2 \
          --frq FRQ \
          --p P \
          --chunksize 50000 \
          --signed-sumstats Z,0 \
          --merge-alleles ${w_hm3} \
          --out ${ID}
    else
      echo "Warning message: please check the input flag for " $ID
    fi
}


ldsc_h2() {
    local ID=$1

    echo "Calculating heritability for ${ID}..."

    ${LDSC}/ldsc.py \
        --h2 ${ID}.sumstats.gz \
        --ref-ld-chr ${ref_ld} \
        --w-ld-chr ${w_ld} \
        --out ${ID}_h2
    
    grep "Total Observed scale h2" ${ID}_h2.log
}

# # PD_2019 ###################################################
# ldsc_munge PD_2019 ${DataPath}/${ID} true Neff
# ldsc_h2 PD_2019 ${DataPath}/${ID}

# # NEA_2021 ###################################################
# ldsc_munge NEA_2021 ${DataPath}/${ID} false 510795
# ldsc_h2 NEA_2021 ${DataPath}/${ID}


# # CEA_2021 ###################################################
# ldsc_munge CEA_2021 ${DataPath}/${ID} false 257700
# ldsc_h2 CEA_2021 ${DataPath}/${ID}

# # CRP_2022 ###################################################
# ldsc_munge CRP_2022 ${DataPath}/${ID} false 575531
# ldsc_h2 CRP_2022 ${DataPath}/${ID}

# # INT_2017 ###################################################
# ldsc_munge INT_2017 ${DataPath}/${ID} false 78308
# ldsc_h2 INT_2017 ${DataPath}/${ID}

# # MDD_2025 ###################################################
# ldsc_munge MDD_2025 ${DataPath}/${ID} true N
# ldsc_h2 MDD_2025 ${DataPath}/${ID}

# # Summary information ########################################
# cd /lustre/home/sww208/GoDMC/HormiseGWAS || exit
# touch collect_h2.log
# for ID in PD_2019 NEA_2021 CEA_2021 CRP_2022 INT_2017 MDD_2025
# do 
#     echo ==================================================== >> collect_h2.log
#     echo $ID >> collect_h2.log
#     grep "Total Observed scale h2" ${DataPath}/${ID}/${ID}_h2.log >> collect_h2.log
#     echo "" >>  collect_h2.log
 
# done


# AD_2021.b37.filtered ###################################################
ldsc_munge AD_2021 true N
ldsc_h2 AD_2021

# T2D_2024.b37.filtered ###################################################
# ldsc_munge T2D_2024 ${DataPath} case N_CAS N_CON
# ldsc_h2 T2D_2024 ${DataPath}

# BMI_2018.b37.filtered ###################################################
ldsc_munge BMI_2018 true N
ldsc_h2 BMI_2018

# BP_2024.b37.filtered ###################################################
${LDSC}/munge_sumstats.py \
  --sumstats BP_2024_MSS.tsv \
  --snp SNP \
  --N-cas-col N_CAS \
  --N-con-col N_CON \
  --a1 A1 \
  --a2 A2 \
  --frq FRQ \
  --p P \
  --chunksize 50000 \
  --signed-sumstats OR,1 \
  --merge-alleles ${w_hm3} \
  --out BP_2024
ldsc_h2 BP_2024

# FG_2021.b37.filtered ###################################################
ldsc_munge FG_2021 true N
ldsc_h2 FG_2021

# FI_2021.b37.filtered ###################################################
ldsc_munge FI_2021 true N
ldsc_h2 FI_2021

# G2H_2021.b37.filtered ###################################################
ldsc_munge G2H_2021 true N
ldsc_h2 G2H_2021

# HBAC1_2021.b37.filtered ###################################################
ldsc_munge HBAC1_2021 true N
ldsc_h2 HBAC1_2021

# HEIGHTb_2022.b37.filtered ###################################################
ldsc_munge HEIGHTb_2022 true N
ldsc_h2 HEIGHTb_2022

# HEIGHTc_2022.b37.filtered ###################################################
ldsc_munge HEIGHTc_2022 true N
ldsc_h2 HEIGHTc_2022

# MD_2025.b37.filtered ###################################################
ldsc_munge MD_2025 case N_CAS NCONTROLS
ldsc_h2 MD_2025

# SCZ_2022.b37.filtered ###################################################
ldsc_munge SCZ_2022 case NCAS NCON
ldsc_h2 SCZ_2022

# T2D_2024_MSS.tsv###################################################
ldsc_munge T2D_2024 case N_CAS N_CON
ldsc_h2 T2D_2024

cd /lustre/home/sww208/GoDMC/HormiseGWAS || exit
touch collect_h2_ann.log
for ID in AD_2021 BP_2024 BMI_2018 FG_2021 FI_2021 G2H_2021 HBAC1_2021 HEIGHTb_2022 HEIGHTc_2022 MD_2025 SCZ_2022 T2D_2024 
do 
    echo ==================================================== >> collect_h2_ann.log
    echo $ID >> collect_h2_ann.log
    grep "Total Observed scale h2" ${DataPath}/${ID}_h2.log >> collect_h2_ann.log
    echo "" >>  collect_h2_ann.log
 
done
