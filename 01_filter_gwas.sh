#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH -A Research_Project-MRC190311 # research project to submit under.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --output=filter.out
#SBATCH --error=filter.err
#SBATCH --job-name=filter
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --ntasks=2
#SBATCH --time=0-10:00:00

############################################################################################
#       Filter GWAS by MAF / INFO
############################################################################################

# for each GWAS, check formats, and write filtered versions to DIRFILT, with suffix indicating build

# use filtering params MAF>1%, INFO>0.6
# use head command to check which columns to specify in filtering

#define the path
RscriptsPath="/lustre/home/sww208/GoDMC/HormiseGWAS/Rscripts"
DataPath="/lustre/home/sww208/GoDMC/GADatasets/OtherGWASstats"
module load R/4.2.1-foss-2022a
# gz function
filter_gz() {
    local ID=$1
    local FILENAME=$2
    local DATAPATH=$3

    cd ${DATAPATH}/${ID} || exit 1

    case "$FILENAME" in
        *.gz) zcat ${FILENAME} | head -n 3 > ${ID}_header.tsv 
              gunzip -c ${FILENAME} > ${ID}.tsv ;;
        *)    head -n 3 ${FILENAME} > ${ID}_header.tsv 
              cp ${FILENAME} ${ID}.tsv ;;
    esac
}

# filter function
filter_gwas(){
    local ID=$1
    local DATAPATH=$2

    cd ${DATAPATH}/${ID} || exit 1

    MAFcol=$(head -n 1 ${ID}_header.tsv | tr '\t' '\n' | nl | grep -E "effect_allele_frequency|MAF|EAF|EAF_HRC" | cut -f 1)
    INFOcol=$(head -n 1 ${ID}_header.tsv | tr '\t' '\n' | nl | grep -E "INFO" | cut -f 1)

    Rscript ${RscriptsPath}/filter_gwas.R \
        "${DATAPATH}/${ID}" \
        "${ID}.tsv" \
        "${ID}" \
        $MAFcol \
        $INFOcol

    #rm "${ID}_header.tsv"
}


# PD_2019 ###################################################
filter_gz PD_2019 GCST009325.tsv ${DataPath}
filter_gwas PD_2019 ${DataPath}

# NEA_2021 ###################################################
filter_gz NEA_2021 GCST90011874_buildGRCh37.tsv.gz ${DataPath}
filter_gwas NEA_2021 ${DataPath}

# CEA_2021 ###################################################
filter_gz CEA_2021 GCST90011875_buildGRCh37.tsv.gz ${DataPath}
filter_gwas CEA_2021 ${DataPath}

# CRP_2022 ###################################################
filter_gz CRP_2022 GCST90029070_buildGRCh37.tsv.gz ${DataPath}
filter_gwas CRP_2022 ${DataPath}

# INT_2017 ###################################################
filter_gz INT_2017 28530673-GCST004364-EFO_0004337-Build37.f.tsv.gz ${DataPath}
filter_gwas INT_2017 ${DataPath}

# MDD_2025 ###################################################
filter_gz MDD_2025 GCST90468123.tsv.gz ${DataPath}
filter_gwas MDD_2025 ${DataPath}

# Sep ###############################################################

# from thessgac
# EA_2016 ###################################################
filter_gz EA_2016 EduYears_Main.txt ${DataPath}
filter_gwas EA_2016 ${DataPath}
awk 'NR==1 {print $0 "\tN"; next} {print $0 "\t293723"}' ${DataPath}/EA_2016/EA_2016.LDSCinput > ${DataPath}/EA_2016/EA_2016.thessgac
awk 'NR==1 {$1="variant_id"; $2="Chromosome"; $3="Position"; $4="effect_allele"; $5="non_effect_allele"; $6="FRQ"; $7="beta"; $8="SE"; $9="Pvalue"; $10="N"} {print}' OFS="\t" ${DataPath}/EA_2016/EA_2016.thessgac > ${DataPath}/EA_2016/EA_2016.LDSCinput

# EA_2022 ###################################################
filter_gz EA_2022 EA4_dominance_excl_23andMe.txt.gz ${DataPath}
mv ${DataPath}/EA_2022/EA_2022.tsv ${DataPath}/EA_2022/EA_2022.tsv.temp
awk '{print $1, $2, $3, $4, $5, $6, $8, $9, $12}' OFS="\t" ${DataPath}/EA_2022/EA_2022.tsv.temp > ${DataPath}/EA_2022/EA_2022.tsv.temp2
awk 'NR==1 {print $0 "\tN"; next} {print $0 "\t765283"}' ${DataPath}/EA_2022/EA_2022.tsv.temp2 > ${DataPath}/EA_2022/EA_2022.tsv
head -3 ${DataPath}/EA_2022/EA_2022.tsv > ${DataPath}/EA_2022/EA_2022_header.tsv
filter_gwas EA_2022 ${DataPath}
mv ${DataPath}/EA_2022/EA_2022.LDSCinput ${DataPath}/EA_2022/EA_2022.thessgac
awk 'NR==1 {$1="variant_id"; $2="Chromosome"; $3="Position"; $4="effect_allele"; $5="non_effect_allele"; $6="FRQ"; $7="beta"; $8="SE"; $9="Pvalue"; $10="N"} {print}' OFS="\t" ${DataPath}/EA_2022/EA_2022.thessgac > ${DataPath}/EA_2022/EA_2022.LDSCinput

# from GIANT
# WC_2015 ###################################################



# WHRBMI_2018 ###################################################
# no extract header

# WHR_2015

# BMI_2018

# BMI_2015

# HEIGHT_2010

# HEIGHT_2022

