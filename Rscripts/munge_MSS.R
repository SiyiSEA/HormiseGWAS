###############################description###################################


#################################packages####################################
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(MungeSumstats))
options(bitmapType = "cairo")


#################################main code####################################
arguments <- commandArgs(T)
setpath <- arguments[1]
filename <- arguments[2]
phenotype <- arguments[3]

setwd(setpath)
print(paste("Working directory set to:", setpath))
GWAresult = read_sumstats(filename)
print(GWAresult[1:3,])

collist = paste(colnames(GWAresult), collapse = " ")


# Check if the data needs compute_z
if ( grepl(c("z|z_score"), collist, ignore.case = TRUE) ) {
  flag_z <- FALSE
  print(paste0(filename, " doesn't need z computation."))
}else{
  flag_z <- TRUE
  print(paste0(filename, " needs z computation."))
}


# Check if the data needs rsid

if ( grepl(c("variant_id|RSID|ID|rsID"), collist, ignore.case = TRUE) ) {
  flag_rsid <- TRUE
  data.table::setnames(GWAresult, c("RSID"), c("variant_id"), skip_absent = T)
  data.table::setnames(GWAresult, c("ID"), c("variant_id"), skip_absent = T)
  print(paste0(filename, " has rsid."))
}else{
  study_list_has_rsid = c("BMI_2018.b37.filtered", "BP_2025.b37.filtered")
  if (filename %in% study_list_has_rsid ){
    flag_rsid <- TRUE
    print(paste0(filename, " has rsid."))
  }else{
    flag_rsid <- FALSE
    print(paste0(filename, " doesn't has rsid.")) 
  }
}


# Check if data needs compute N
if ( grepl(c("n|N_cases|N_controls|Nca|Nco|sample_size|ncontrols|ncases|NCAS|NCON"), collist, ignore.case = TRUE) ) {
  flag_n <- 0
  data.table::setnames(GWAresult, c("sample_size"), c("N"), skip_absent = T)
  print(paste0(filename, " doesn't need n computation."))
}else{
  flag_n <- "ldsc"
  print(paste0(filename, " needs n computation."))
}
if (filename %in% c("T2D_2024.b37.filtered.SNP", "BP_2025.b37.filtered")){
  flag_n <- "ldsc"
  print(paste0("Actually!! ",filename, " needs n computation."))
}

# set A1 and A2 in a proper name
if ("effect_allele" %in% colnames(GWAresult)){
  print(paste0(filename, " has effect_allele."))
}else{
  print(paste0("Setting effect_allele for ", filename))
  if (grepl(c("A1|tested_allele|EffectAllele|Effect_allele"), collist, ignore.case = TRUE)){
    data.table::setnames(GWAresult, c("A1"), c("effect_allele"), skip_absent = T)
    data.table::setnames(GWAresult, c("tested_allele"), c("effect_allele"), skip_absent = T)
    data.table::setnames(GWAresult, c("Tested_Allele"), c("effect_allele"), skip_absent = T)
    data.table::setnames(GWAresult, c("EffectAllele"), c("effect_allele"), skip_absent = T)
    data.table::setnames(GWAresult, c("EFFECT_ALLELE"), c("effect_allele"), skip_absent = T)
    data.table::setnames(GWAresult, c("Effect_allele"), c("effect_allele"), skip_absent = T)
  }else{
    message("Warning: there is no colname for the effect allele, please check the data!")
  }
}

if ("non_effect_allele" %in% colnames(GWAresult)){
  print(paste0(filename, " has non_effect_allele."))
}else{
  print(paste0("Setting non_effect_allele for ", filename ))
  if (grepl(c("other_allele|A2|NonEffectAllele"), collist, ignore.case = TRUE)){
    data.table::setnames(GWAresult, c("other_allele"), c("non_effect_allele"), skip_absent = T)
    data.table::setnames(GWAresult, c("Other_Allele"), c("non_effect_allele"), skip_absent = T)
    data.table::setnames(GWAresult, c("OTHER_ALLELE"), c("non_effect_allele"), skip_absent = T)
    data.table::setnames(GWAresult, c("A2"), c("non_effect_allele"), skip_absent = T)
    data.table::setnames(GWAresult, c("NonEffectAllele"), c("non_effect_allele"), skip_absent = T)
    data.table::setnames(GWAresult, c("Other_allele"), c("non_effect_allele"), skip_absent = T)
  }else{
    message("Warning: there is no colname for the non_effect_allele, please check the data!")
  }
}


# set effect allele frequency in a proper name
if ("FRQ" %in% colnames(GWAresult)){
  print(paste0(filename, " has FRQ."))
}else{
  print(paste0("Setting FRQ for ", filename ))
  if (grepl(c("effect_allele_frequency|MAF|EAF|Freq_Tested_Allele_in_HRS|HRC_FRQ_A1|EFFECT_ALLELE_FREQ"), collist, ignore.case = TRUE)){
    data.table::setnames(GWAresult, c("effect_allele_frequency"), c("FRQ"), skip_absent = T)
    data.table::setnames(GWAresult, c("MAF"), c("FRQ"), skip_absent = T)
    data.table::setnames(GWAresult, c("EAF"), c("FRQ"), skip_absent = T)
    data.table::setnames(GWAresult, c("HRC_FRQ_A1"), c("FRQ"), skip_absent = T)
    data.table::setnames(GWAresult, c("EFFECT_ALLELE_FREQ"), c("FRQ"), skip_absent = T)
    data.table::setnames(GWAresult, c("Freq_Tested_Allele_in_HRS"), c("FRQ"), skip_absent = T)
  }else{
    message("Warning: there is no colname for the frq, please check the data!")
  }
}

# set BETA
if ("est" %in% colnames(GWAresult)){
  data.table::setnames(GWAresult, c("est"), c("beta"))
}


print(GWAresult[1:3,])

reformatted <- MungeSumstats::format_sumstats(path = GWAresult,
                                          ref_genome="GRCH37",
                                          on_ref_genome=TRUE,
                                          compute_z=flag_z, 
                                          compute_n=flag_n,
                                          snp_ids_are_rs_ids = flag_rsid,
                                          bi_allelic_filter = FALSE,
                                          flip_frq_as_biallelic = TRUE,
                                          drop_indels = TRUE,
                                          nThread = 6,
                                          save_path = paste0(setpath, "/",phenotype,"_MSS.tsv"),
                                          log_mungesumstats_msgs = TRUE,
                                          force_new=TRUE,
                                          save_format='LDSC',
                                          log_folder = setpath)

print(head(reformatted))
save(flag_z, flag_n, flag_rsid, reformatted, file = paste0(phenotype,"_MSS.RData"))