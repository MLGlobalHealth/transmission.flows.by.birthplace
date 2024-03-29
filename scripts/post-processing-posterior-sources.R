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

infile.seq <-	file.path(args_dir$indir, 'Data', 'data_220331/SHM_2201_ROADMAP_220331_tblLAB_seq.rda')
infile.meta <- file.path(args_dir$indir, args_dir$analysis, 'misc', '220713_sequence_labels.rda')
infile.po.tpairprob <- paste0(outfile.base,'-stanout-tpairprobw-gqs.RDS')
infile.sampling.prob <- paste0(outfile.base,'-sampling_prob_byyear_cases','.RDS')

cat('\nReading patient metadata...')
load(infile.meta)
load(infile.seq)

cat('\nReading sampling probabilities...')
spy <- readRDS(file=infile.sampling.prob)
cat('\nReading posterior transmission pair probabilities...')
tprob <- readRDS(file=infile.po.tpairprob)
#samples <- rstan::extract(fit, inc_warmup = FALSE)
tprob <- data.table(reshape::melt(tprob$tpair_prob_w))

cat(" \n --------------------------------  add birthplace data to pairs -------------------------------- \n")

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
# estimate sources by birthplace ----
cat(" \n --------------------------------  estimate sources by birthplace -------------------------------- \n")

po <- data.table(tprob)
#setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
#po <- melt(po, id.vars = c('chain','iteration','draw'))
#po <- data.table(po)
setnames(po,c('iterations','Var.2'),c('draw','PAIR_ID'))
#po[, PAIR_ID := as.integer(gsub(paste0('tpair_prob_w\\[([0-9]+)\\]'),'\\1',as.character(variable)))]
tmp <- subset(do, select = c('PAIR_ID','FROM_BPLACE','TO_BPLACE','YEAR_OF_INF_EST'))
po <- merge(po, tmp, by = 'PAIR_ID')
po <- merge(po, subset(spy,select=c('LOC_BIRTH_POS','YEAR_OF_INF_EST','psi')),
            by.x=c('TO_BPLACE','YEAR_OF_INF_EST'),by.y=c('LOC_BIRTH_POS','YEAR_OF_INF_EST'))
po <- po[, list(value = sum(value/psi)), by = c('draw','FROM_BPLACE')]
tmp <- po[, list(total = sum(value)), by = c('draw')]
po <- merge(po, tmp, by = 'draw')
po[, paf := value/total]
saveRDS(po,file=paste0(outfile.base,'-flows_frombplace_adjusted_samplingbias_mcsamples','.RDS'))


# stratified flows by birthplace of case and source ----
cat(" \n --------------------------------  plot adjusted flows by birthplace -------------------------------- \n")
## % flows from each region PER region of birth of cases (i.e. sum to one per regino of birth of recipient)

po <- data.table(tprob)
setnames(po,c('iterations','Var.2'),c('draw','PAIR_ID'))
tmp <- subset(do, select = c('PAIR_ID','FROM_BPLACE','TO_BPLACE','YEAR_OF_INF_EST'))
po <- merge(po, tmp, by = 'PAIR_ID')
po <- merge(po, subset(spy,select=c('LOC_BIRTH_POS','YEAR_OF_INF_EST','psi')),
            by.x=c('TO_BPLACE','YEAR_OF_INF_EST'),by.y=c('LOC_BIRTH_POS','YEAR_OF_INF_EST'))
po <- po[, list(value = sum(value/psi)), by = c('draw','FROM_BPLACE','TO_BPLACE')]
tmp <- po[, list(total = sum(value)), by = c('draw','TO_BPLACE')]
po <- merge(po, tmp, by = c('draw','TO_BPLACE'))
po[, paf := value/total]

saveRDS(po,file=paste0(outfile.base,'-stratified_flows_frombplace_adjusted_samplingbias_mcsamples','.RDS'))

## get flows from group a to group b out of total flows ----
cat(" \n --------------------------------  plot flows from group a to group b out of total flows -------------------------------- \n")

po <- data.table(tprob)
setnames(po,c('iterations','Var.2'),c('draw','PAIR_ID'))
tmp <- subset(do, select = c('PAIR_ID','FROM_BPLACE','TO_BPLACE','YEAR_OF_INF_EST'))
po <- merge(po, tmp, by = 'PAIR_ID')
po <- merge(po, subset(spy,select=c('LOC_BIRTH_POS','YEAR_OF_INF_EST','psi')),
            by.x=c('TO_BPLACE','YEAR_OF_INF_EST'),by.y=c('LOC_BIRTH_POS','YEAR_OF_INF_EST'))
po <- po[, list(value = sum(value/psi)), by = c('draw','FROM_BPLACE','TO_BPLACE')]
tmp <- po[, list(total = sum(value)), by = c('draw')]
po <- merge(po, tmp, by = c('draw'))
po[, paf := value/total]

saveRDS(po,file=paste0(outfile.base,'-adjusted_flows_atob_samplingofcases_bplacecase_bplacesrc','.RDS'))

do[, FROM_BPLACE:= factor(FROM_BPLACE,
                          levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\n Caribbean','E. & C. Europe','MENA','Other'),
                          labels=c('Netherlands','W.Europe,\nN.America & Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\nNon-Dutch Caribbean','E. & C. Europe','MENA','Other'))]
do[, TO_BPLACE:= factor(TO_BPLACE,
                        levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                 'S. America &\n Caribbean','E. & C. Europe','MENA','Other'),
                        labels=c('Netherlands','W.Europe,\nN.America & Oceania','Suriname &\nDutch Caribbean',
                                 'S. America &\nNon-Dutch Caribbean','E. & C. Europe','MENA','Other'))]

breaks_to <- do[, list(N_TO=length(unique(TO_SEQUENCE_ID))),by=c('TO_BPLACE')]
breaks_from <- do[, list(N_FROM=length(unique(FROM_SEQUENCE_ID))),by=c('FROM_BPLACE')]
breaks_to <- breaks_to[order(TO_BPLACE),]
breaks_from <- breaks_from[order(FROM_BPLACE),]

breaks_to[, pos_to:= cumsum(N_TO) - N_TO/2 ]
breaks_from[, pos_from:= cumsum(N_FROM) - N_FROM/2]

saveRDS(breaks_from,file=paste0(outfile.base,'-flows_breaks_from','.RDS'))
saveRDS(breaks_to,file=paste0(outfile.base,'-flows_breaks_to','.RDS'))
