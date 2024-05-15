## preamble ----
require(data.table)

if (0)
{
  args_dir <- list(
    source_dir = '~/Documents/GitHub/transmission.flows.by.birthplace',
    #indir = '~/Box\ Sync/Roadmap/source_attribution',
    indir = '~/Box\ Sync/Roadmap',
    analysis = 'analysis_220713',
    pairs.dir = '/Users/alexb/Documents/GitHub/transmission.flows.by.birthplace/out_Amsterdam/published_data_MSM-2010_2021',
    outdir = '/Users/alexb/Documents/GitHub/transmission.flows.by.birthplace/out_Amsterdam/mm_bgUnif_piGP_221027b-published_data_MSM-2010_2021-860418',
    clock_model = '/Users/alexb/Box Sync/Roadmap/source_attribution/molecular_clock/hierarchical',
    stanModelFile = 'mm_bgUnif_piGP_221027b',
    job_tag = 'published_data_MSM-2010_2021',
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

infile.sampling.prob <- file.path(args_dir$source_dir,'data','sampling_prob_byyear_bplace.RDS')

tmp <- paste0(outfile.base,'-fitted_stan_model.rds')
cat("\n Read fitted dat ", tmp , "\n")
model_fit <- readRDS(file = tmp)

cat('\nReading sampling probabilities...')
spy <- readRDS(file=infile.sampling.prob)

# summarise flows by two-year intervals ----

cat('\nReading posterior transmission pair probabilities...')
po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = 'tpair_prob_w' #tpair_prob_w
)
po <- data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po <- data.table(po)
po[, PAIR_ID := as.integer(gsub(paste0('tpair_prob_w\\[([0-9]+)\\]'),'\\1',as.character(variable)))]
tmp <- subset(do, select = c('PAIR_ID','FROM_BPLACE','TO_BPLACE','FROM_SEQUENCE_ID','TO_SEQUENCE_ID','FROM_MIGRANT','TO_MIGRANT','YEAR_OF_INF_EST'))
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
