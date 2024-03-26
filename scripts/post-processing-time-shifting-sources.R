## preamble ----
require(data.table)

if (0)
{
  args_dir <- list(
    source_dir = '~/Documents/GitHub/transmission.flows.by.birthplace',
    #indir = '~/Box\ Sync/Roadmap/source_attribution',
    indir = '~/Box\ Sync/Roadmap',
    analysis = 'analysis_220713',
    pairs.dir = '/Users/alexb/Documents/GitHub/transmission.flows.by.birthplace/out_Amsterdam/agegps_updated_criteria_MSM-2010_2022',
    outdir = '/Users/alexb/Documents/GitHub/transmission.flows.by.birthplace/out_Amsterdam/mm_bgUnif_piGP_221027b-agegps_sensanalysis_210216_MSM-618873',
    clock_model = '/Users/alexb/Box Sync/Roadmap/source_attribution/molecular_clock/hierarchical',
    stanModelFile = 'mm_bgUnif_piGP_221027b',
    scenario = 15,
    reps = 1,
    rep = 1,
    simulate_data = T,
    job_tag = 'agegps_sensanalysis_210216_MSM-2010_2022',
    undiagnosed = '/Users/alexb/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Roadmap/sources/ethnicity_analysis/undiagnosed/undiagnosed_211102-cohort_2010_2015',
    overwrite = 1
  )
}

## command line parsing if any
args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0)
{
  stopifnot(args_line[[1]]=='-source_dir')
  stopifnot(args_line[[3]]=='-indir')
  stopifnot(args_line[[5]]=='-outdir')
  stopifnot(args_line[[7]]=='-stanModelFile')
  stopifnot(args_line[[9]]=='-job_tag')
  stopifnot(args_line[[11]]=='-analysis')

  args_dir <- list()
  args_dir[['source_dir']] <- args_line[[2]]
  args_dir[['indir']] <- args_line[[4]]
  args_dir[['outdir']] <- args_line[[6]]
  args_dir[['stanModelFile']] <- args_line[[8]]
  args_dir[['job_tag']] <- args_line[[10]]
  args_dir[['analysis']] <- args_line[[12]]
}
args_dir

cat(" \n --------------------------------  load data -------------------------------- \n")

## read stanin
cat('\nReading Stan input data...')
infile.stanin <- list.files(args_dir$outdir, pattern=paste0('_stanin.RData$'), recursive=TRUE)[1]
stopifnot(length(infile.stanin)>0)
stopifnot(length(infile.stanin)<=1)
tmp <- load(file.path(args_dir$outdir, infile.stanin))
stopifnot(c('args','stan_data')%in%tmp)

outfile.base <- paste0(args_dir$outdir, "/",
                       args_dir$stanModelFile , "-", args_dir$job_tag)


infile.meta <- file.path(args_dir$indir, args_dir$analysis, 'misc', '220713_sequence_labels.rda')
infile.po.tpairprob <- paste0(outfile.base,'-stanout-tpairprobw-gqs.RDS')
infile.sampling.prob <- paste0(outfile.base,'-sampling_prob_byyear_cases','.RDS')

cat('\nReading sampling probabilities...')
spy <- readRDS(file=infile.sampling.prob)
cat('\nReading posterior transmission pair probabilities...')
tprob <- readRDS(file=infile.po.tpairprob)
#samples <- rstan::extract(fit, inc_warmup = FALSE)
tprob <- data.table(reshape::melt(tprob$tpair_prob_w))

## load ethnicity data ----

cat('\nReading patient metadata...')
load(infile.meta)

dind <- data.table(dind)

do <- unique(do)
do <- merge(do,subset(dind,select=c('PATIENT','LOC_BIRTH')),by.x='FROM_SEQUENCE_ID',by.y='PATIENT',all.x=T)
setnames(do,'LOC_BIRTH','FROM_LOC_BIRTH')
do <- merge(do,subset(dind,select=c('PATIENT','LOC_BIRTH')),by.x='TO_SEQUENCE_ID',by.y='PATIENT',all.x=T)
setnames(do,'LOC_BIRTH','TO_LOC_BIRTH')

do <- unique(do)

do[, FROM_BPLACE:="Other"]
do[FROM_LOC_BIRTH %in% c("WEurope","NorthAm","Oceania"), FROM_BPLACE:="W.Europe,\nN.America,Oceania"]
do[FROM_LOC_BIRTH %in% c("EEurope", "CEurope"), FROM_BPLACE:="E. & C. Europe"]
do[FROM_LOC_BIRTH %in% c("LaAmCar"), FROM_BPLACE:="S. America &\n Caribbean"]
do[FROM_LOC_BIRTH %in% c("DutchCarSuriname"), FROM_BPLACE:="Suriname &\nDutch Caribbean"]
do[FROM_LOC_BIRTH %in% c("MENA"), FROM_BPLACE:="MENA"]
do[FROM_ORIGIN=="NL", FROM_BPLACE:="Netherlands"]

do[, TO_BPLACE:="Other"]
do[TO_LOC_BIRTH %in% c("WEurope","NorthAm","Oceania"), TO_BPLACE:="W.Europe,\nN.America,Oceania"]
do[TO_LOC_BIRTH %in% c("EEurope", "CEurope"), TO_BPLACE:="E. & C. Europe"]
do[TO_LOC_BIRTH %in% c("LaAmCar"), TO_BPLACE:="S. America &\n Caribbean"]
do[TO_LOC_BIRTH %in% c("DutchCarSuriname"), TO_BPLACE:="Suriname &\nDutch Caribbean"]
do[TO_LOC_BIRTH %in% c("MENA"), TO_BPLACE:="MENA"]
do[TO_ORIGIN=="NL", TO_BPLACE:="Netherlands"]

do[, FROM_BPLACE:= factor(FROM_BPLACE,
                          levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\n Caribbean','E. & C. Europe','MENA','Other'))]
do[, TO_BPLACE:= factor(TO_BPLACE,
                        levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                 'S. America &\n Caribbean','E. & C. Europe','MENA','Other'))]

# summarise flows by two-year intervals ----

po <- data.table(tprob)
setnames(po,c('iterations','Var.2'),c('draw','PAIR_ID'))
tmp <- subset(do, select = c('PAIR_ID','FROM_BPLACE','TO_BPLACE','FROM_COUNTRY','TO_COUNTRY','FROM_SEQUENCE_ID','TO_SEQUENCE_ID','FROM_MIGRANT','TO_MIGRANT','YEAR_OF_INF_EST'))
po <- merge(po, tmp, by = 'PAIR_ID')
po <- merge(po, subset(spy,select=c('LOC_BIRTH_POS','YEAR_OF_INF_EST','psi')),
            by.x=c('TO_BPLACE','YEAR_OF_INF_EST'),by.y=c('LOC_BIRTH_POS','YEAR_OF_INF_EST'))
po[, YEAR_GP:= cut(YEAR_OF_INF_EST, breaks=c(2010,2012,2014,2016,2018,2020,2022),
                   labels=c('2010-2011','2012-2013','2014-2015','2016-2017','2018-2019','2020-2021'),include.lowest=T,right=F)]
# print sample sizes
cat('\nNumber of incident cases per 2-year interval:\n')
po[, list(N_cases=length(unique(TO_SEQUENCE_ID))),by='YEAR_GP']
cat('\nNumber of pairs per 2-year interval:\n')
po[, list(N_pairs=length(unique(PAIR_ID))),by='YEAR_GP']

po <- po[, list(value = sum(value/psi)), by = c('draw','FROM_MIGRANT','YEAR_GP')]
tmp <- po[, list(total = sum(value)), by = c('draw','YEAR_GP')]
po <- merge(po, tmp, by = c('draw','YEAR_GP'))
po[, paf := value/total]

saveRDS(po,file=paste0(outfile.base,'-adjusted_flows_mwmb_by2years_samplingofcases_mcsamples','.RDS'))
