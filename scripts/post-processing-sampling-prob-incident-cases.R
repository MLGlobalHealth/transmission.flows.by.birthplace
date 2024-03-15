## preamble ----
require(data.table)
require(ggplot2)
require(ggsci)
require(scales)
require(grid)
require(dplyr)
require(tidyr)
require(lubridate)

if (1)
{
  args <- list(
    source_dir = '~/Documents/GitHub/transmission.flows.by.birthplace',
    #indir = '~/Box\ Sync/Roadmap/source_attribution',
    indir = '~/Box\ Sync/Roadmap',
    analysis = 'analysis_220713',
    pairs.dir = '/Users/alexb/Documents/GitHub/source.attr.with.infection.time.fork/out_Amsterdam/agegps_updated_criteria_MSM-2010_2022',
    #outdir = '/Users/alexb/Documents/GitHub/source.attr.with.infection.time.fork/out_Amsterdam/mm_bgUnif_piGP_221027b-agegps_TE16_MSM-2010_2022-1665619',
    outdir = '/Users/alexb/Documents/GitHub/source.attr.with.infection.time.fork/out_Amsterdam/mm_bgUnif_piGP_221027b-agegps_sensanalysis_210216_MSM-618873',
    clock_model = '/Users/alexb/Box Sync/Roadmap/source_attribution/molecular_clock/hierarchical',
    stanModelFile = 'mm_bgUnif_piGP_221027b',
    scenario = 15,
    reps = 1,
    rep = 1,
    simulate_data = T,
    job_tag = 'agegps_sensanalysis_210216_MSM-2010_2022',
    undiagnosed = '/Users/alexb/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Roadmap/sources/ethnicity_analysis/undiagnosed/undiagnosed_211102-cohort_2010_2015',
    job_tag_undiag = 'cohort_2010_2015'
  )
}

## command line parsing if any
args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0)
{
  stopifnot(args_line[[1]]=='-source_dir')
  stopifnot(args_line[[3]]=='-stanModelFile')
  stopifnot(args_line[[5]]=='-indir')
  stopifnot(args_line[[7]]=='-outdir')
  stopifnot(args_line[[9]]=='-job_tag')
  stopifnot(args_line[[11]]=='-undiagnosed')
  stopifnot(args_line[[13]]=='-job_tag_undiag')
  stopifnot(args_line[[15]]=='-analysis')

  args_dir <- list()
  args_dir[['source_dir']] <- args_line[[2]]
  args_dir[['stanModelFile']] <- args_line[[4]]
  args_dir[['indir']] <- args_line[[6]]
  args_dir[['outdir']] <- args_line[[8]]
  args_dir[['job_tag']] <- args_line[[10]]
  args_dir[['undiagnosed']] <- args_line[[12]]
  args_dir[['job_tag_undiag']] <- args_line[[14]]
  args_dir[['analysis']] <- args_line[[16]]
}
args

source(file.path(args_dir$source_dir, 'R', 'functions.R'))

cat(" \n --------------------------------  with arguments -------------------------------- \n")

infile.seq <-	file.path(args_dir$indir, 'Data', 'data_220331/SHM_2201_ROADMAP_220331_tblLAB_seq.rda')
infile.bas <- file.path(args_dir$indir, 'Data', 'data_220331','SHM_2201_ROADMAP_220331_tblBAS.csv')
infile.meta <- file.path(args_dir$indir, args_dir$analysis, 'misc', '220713_sequence_labels.rda')

## read stanin
cat('\nReading Stan input data...')
infile.stanin <- list.files(args_dir$outdir, pattern=paste0('_stanin.RData$'), recursive=TRUE)[1]
stopifnot(length(infile.stanin)>0)
stopifnot(length(infile.stanin)<=1)
tmp <- load(file.path(args_dir$outdir, infile.stanin))
stopifnot(c('args','stan_data')%in%tmp)

# calculate sampling prop for incident cases ----

## load infection date data ----
cat('\nReading infection date estimates...')
dinf <- data.table(read.csv(file.path(args_dir$indir,'transmission_sources','Infection_date_est_rec.csv')))
setnames(dinf,c("id",'estsctodiagMedian','estsctodiagLL','estsctodiagUL'),c("TO_SEQUENCE_ID",'SER_TO_DIAG','SER_TO_DIAG_LL','SER_TO_DIAG_UL'))
dinf <- unique(dinf)
dinf[,DIAGNOSIS_DATE:= as.Date(dinf[,hiv_pos_d],format="%Y-%m-%d")]
dinf[,DIAGNOSIS_DATE_N:= hivc.db.Date2numeric(dinf[,DIAGNOSIS_DATE])]
dinf[,EST_INF_DATE:= DIAGNOSIS_DATE_N-SER_TO_DIAG]
dinf[,EST_INF_DATE:= format(date_decimal(EST_INF_DATE), "%Y-%m-%d")]
dinf[,YEAR_OF_INF_EST := year(EST_INF_DATE)]

### merge in patient metadata ----
cat('\nReading patient metadata...')
load(infile.meta)
load(infile.seq)

dind <- data.table(dind)
dind[, SEQ:= PATIENT %in% ds$PATIENT]
dinf <- merge(dinf,subset(dind,select=c('PATIENT','CITY','TRANSM','LOC_BIRTH','ORIGIN','SEQ')),
              by.x='TO_SEQUENCE_ID',by.y='PATIENT',all.x=T)
dbas <- fread(infile.bas)
dbas[, DEATH_D := as.Date(DEATH_D,format="%Y-%m-%d")]
dinf <- merge(dinf,subset(dbas,select=c('PATIENT','DEATH_D')),by.x='TO_SEQUENCE_ID',by.y='PATIENT',all.x=T)

## add world regions ----
cat('\nAdding world regions...')

dinf[, LOC_BIRTH_POS:="Other"]
dinf[LOC_BIRTH %in% c("WEurope","NorthAm","Oceania"), LOC_BIRTH_POS:="W.Europe,\nN.America,Oceania"]
dinf[LOC_BIRTH %in% c("EEurope", "CEurope"), LOC_BIRTH_POS:="E. & C. Europe"]
dinf[LOC_BIRTH %in% c("LaAmCar"), LOC_BIRTH_POS:="S. America &\n Caribbean"]
dinf[LOC_BIRTH %in% c("DutchCarSuriname"), LOC_BIRTH_POS:="Suriname &\nDutch Caribbean"]
dinf[LOC_BIRTH %in% c("MENA"), LOC_BIRTH_POS:="MENA"]
dinf[ORIGIN=="NL", LOC_BIRTH_POS:="Netherlands"]

# calculate number diagnosed from each birth region and number with a sequence

cat('\nLoad posterior from undiagnosed model...')

samples <- readRDS(file=file.path(args_dir$undiagnosed, paste0('samples_',args_dir$job_tag_undiag,"_","MSM",'.rds')))

cat('\nLoading geographic region mapping...')

dmap <- readRDS(file=file.path(args_dir$undiagnosed, paste0("mapping_georeg_id.RDS")))

shape_msm <- data.table(reshape::melt(samples$wb_shape_grp))
setnames(shape_msm,c('iterations','Var.2'),c('iter','mg'))
shape_msm[, trsm:='MSM']
shape_msm[, par:='shape']

scale_msm <- data.table(reshape::melt(samples$wb_scale_grp))
setnames(scale_msm,c('iterations','Var.2'),c('iter','mg'))
scale_msm[, trsm:='MSM']
scale_msm[, par:='scale']

ds <- rbind(shape_msm,scale_msm)
ds <- dcast(ds,trsm+mg+iter~par,value.var="value")

dat <- tidyr::crossing(year=seq(1980,2021,1),month=seq(1,12,1))
#dat <- tidyr::crossing(year=seq(2010,2021,1))
ds <- merge(dat,ds,all=T)
ds <- data.table(ds)
ds[, time:=(2022+(1/12))-(year + (month/12))] # +(1/12) to add 1 month so we count until end Dec 2018/start of Jan 2019. but I think it makes sense someone infected Dec 2018 has 0 prob of being diagnosed by end of 2018 due to time taken to detect virus
ds[, p:=1 - pweibull(time,shape=shape,scale=scale)]

mean_y <- ds[, list(p=quantile(p,prob=c(0.025,0.5,0.975)),
                    qlabel=c('p0.025','p0.5','p0.975')),
             by=c('trsm','mg','year')] # summarise quantiles for each year
mean_y <- merge(mean_y,dmap,by.x='mg',by.y='mgid')
mean_y <- dcast(mean_y,trsm+mwmb+year~qlabel,value.var=c("p"))
saveRDS(mean_y,file=paste0(outfile.base,'-mean_undiagnosed_byyear_sources','.RDS'))

## calculate sp by year for cases ----
dat <- tidyr::crossing(year=seq(2010,2021,1),month=seq(1,12,1))
#dat <- tidyr::crossing(year=seq(2010,2021,1))
ds <- merge(dat,ds,all=T)
ds <- data.table(ds)
ds[, time:=(2022+(1/12))-(year + (month/12))] # +(1/12) to add 1 month so we count until end Dec 2018/start of Jan 2019. but I think it makes sense someone infected Dec 2018 has 0 prob of being diagnosed by end of 2018 due to time taken to detect virus
ds[, p:=1 - pweibull(time,shape=shape,scale=scale)]

mean_y <- ds[, list(p=quantile(p,prob=c(0.025,0.5,0.975)),
                    qlabel=c('p0.025','p0.5','p0.975')),
             by=c('trsm','mg','year')] # summarise quantiles for each year
mean_y <- merge(mean_y,dmap,by.x='mg',by.y='mgid')
mean_y <- dcast(mean_y,trsm+mwmb+year~qlabel,value.var=c("p"))


spy <- dinf[YEAR_OF_INF_EST>=2010, list(N=length(TO_SEQUENCE_ID),
                                        N_seq = length(TO_SEQUENCE_ID[SEQ==T])),
            by=c('LOC_BIRTH_POS','YEAR_OF_INF_EST')]
# add rows for years with missing values for each birth region? would need to make the N=0.0001 or similar..
spy <- merge(spy,mean_y,by.x=c('LOC_BIRTH_POS','YEAR_OF_INF_EST'),by.y=c('mwmb','year'),all.x=T)
spy[, N_inf:= N/(1-p0.5)]
spy[, N_inf_CL:=  N/(1-p0.025)]
spy[, N_inf_UL:=  N/(1-p0.975)]

#spy <- spy[, list(N=sum(N),N_inf=sum(N_inf),N_seq=sum(N_seq)),by=c('LOC_BIRTH_POS')]
# calculate sampling prob among incident cases
spy[, N_seq:= as.numeric(N_seq)]
spy[N_seq==0, N_seq:= 0.1]
spy[, psi:= N_seq/N_inf]
spy[, psi_CL:=  N_seq/N_inf_CL]
spy[, psi_CU:=  N_seq/N_inf_UL]

spy[, LOC_BIRTH_POS:= factor(LOC_BIRTH_POS,
                             levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                      'S. America &\n Caribbean','E. & C. Europe','MENA','Other'))]
saveRDS(spy,file=paste0(outfile.base,'-sampling_prob_byyear_cases','.RDS'))
