
## preamble ----
require(data.table)  # data mangling
require(ggplot2)
require(ggsci)
require(scales)
require(grid)
require(ggpubr)
require(ggExtra)
require(cowplot)
require(Hmisc)
library(networkD3)
library(htmlwidgets)
require(dplyr)
require(tidyr)
source('R/functions.R')

if (1)
{
  args <- list(
    source_dir = '~/Documents/GitHub/source.attr.with.infection.time.fork',
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
    job_tag = 'agegps_sensanalysis_210216_MSM-2010_2022'
    undiagnosed = '/Users/alexb/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Roadmap/sources/ethnicity_analysis/undiagnosed/undiagnosed_211102-cohort_2010_2015'
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
  stopifnot(args_line[[11]]=='-scenario')
  stopifnot(args_line[[13]]=='-rep')
  stopifnot(args_line[[15]]=='-weights')

  args <- list()
  args[['source_dir']] <- args_line[[2]]
  args[['stanModelFile']] <- args_line[[4]]
  args[['indir']] <- args_line[[6]]
  args[['outdir']] <- args_line[[8]]
  args[['job_tag']] <- args_line[[10]]
  args[['scenario']] <- as.integer(args_line[[12]])
  args[['rep']] <- as.integer(args_line[[14]])
  args[['weights']] <- as.integer(args_line[[16]])
}
args

args$analysis = 'analysis_220713'
args$indir = '~/Box\ Sync/Roadmap'

cat(" \n --------------------------------  with arguments -------------------------------- \n")

infile.seq <-	file.path(args$indir, 'Data', 'data_220331/SHM_2201_ROADMAP_220331_tblLAB_seq.rda')
infile.bas <- file.path(args$indir, 'Data', 'data_220331','SHM_2201_ROADMAP_220331_tblBAS.csv')
infile.meta <- file.path(args$indir, args$analysis, 'misc', '220713_sequence_labels.rda')


## load incidence cohort ----
cat(" \n --------------------------------  load incidence cohort -------------------------------- \n")

load(infile.seq)
load(infile.meta)
dind <- data.table(dind)
dind[, SEQ:= PATIENT %in% ds$PATIENT]
dind <- unique(dind)

## load infection time estimates and metadata ----

dinf <- data.table(read.csv(file.path('data_Ams',analysis,'Infection_date_est_rec.csv')))
setnames(dinf,c("id",'estsctodiagMedian','estsctodiagLL','estsctodiagUL'),c("SEQUENCE_ID",'SER_TO_DIAG','SER_TO_DIAG_LL','SER_TO_DIAG_UL'))
dinf <- unique(dinf)
dinf <- merge(dinf,subset(dind,select=c('PATIENT','CITY','SEQ','TRANSM','BIRTH_CNTRY','LOC_BIRTH')),
              by.x='SEQUENCE_ID',by.y='PATIENT',all.x=T)
dinf[, SEQ:= SEQUENCE_ID %in% ds$PATIENT]

# calculate infection date
dinf[,DIAGNOSIS_DATE:= as.Date(dinf[,hiv_pos_d],format="%Y-%m-%d")]
dinf[,DIAGNOSIS_DATE_N:= hivc.db.Date2numeric(dinf[,DIAGNOSIS_DATE])]
dinf[,EST_INF_DATE:= DIAGNOSIS_DATE_N-SER_TO_DIAG]
dinf[,EST_INF_DATE:= format(date_decimal(EST_INF_DATE), "%Y-%m-%d")]
dinf[,YEAR_OF_INF_EST := year(EST_INF_DATE)]

dinf <- unique(subset(dinf,CITY=='Amsterdam' & TRANSM=='MSM' & YEAR_OF_INF_EST >= 2010 & YEAR_OF_INF_EST<2022))
dinf <- merge(dinf,subset(dbas,select=c('PATIENT','BIRTH_D')),by.x='SEQUENCE_ID',by.y='PATIENT',all.x=T)
dinf$EST_INF_DATE <- as.Date(ISOdate(dinf$EST_INF_DATE, 1, 1))
dinf$BIRTH_D <- as.Date(ISOdate(dinf$BIRTH_D, 1, 1))
dinf[,TO_AGE:= as.numeric(EST_INF_DATE-BIRTH_D)/365]

# reclassify birth regions
dinf[, LOC_BIRTH_POS:="Other"]
dinf[LOC_BIRTH %in% c("WEurope","NorthAm","Oceania"), LOC_BIRTH_POS:="W.Europe,\nN.America,Oceania"]
dinf[LOC_BIRTH %in% c("EEurope", "CEurope"), LOC_BIRTH_POS:="E. & C. Europe"]
dinf[LOC_BIRTH %in% c("LaAmCar"), LOC_BIRTH_POS:="S. America &\n Caribbean"]
dinf[LOC_BIRTH %in% c("DutchCarSuriname"), LOC_BIRTH_POS:="Suriname &\nDutch Caribbean"]
dinf[LOC_BIRTH %in% c("MENA"), LOC_BIRTH_POS:="MENA"]
dinf[BIRTH_CNTRY=="Netherlands", LOC_BIRTH_POS:="Netherlands"]

# count number of cases
sp <- dinf[YEAR_OF_INF_EST>=2010, list(N=length(SEQUENCE_ID)),
          by='LOC_BIRTH_POS']

## load estimates of proportion undiagnosed ----

ds <- readRDS(file=file.path(args$undiagnosed,paste0('p_undiagnosed_average_cohort_2010_2015_MSM.RDS')))
sp <- merge(sp,ds,by.x='LOC_BIRTH_POS',by.y='migrant_group')

## calculate total infected
sp[, N_inf:= round(N/(1-p0.5))]

## load VL data and estimate unsuppressed ----
cat(" \n --------------------------------  get unsuppressed by year -------------------------------- \n")
# should be out of OBSERVED COHORT (not undiagnosed)
## OR: should we include the undiagnosed in the numerator too, since presumably they are unsuppressed?
# probably doesn't matter because if we add them to numerator and denominator we add no information, but we may have less uncertainty (bad)

infile.rna <-	file.path(indir_data, 'Data', 'data_220331/SHM_2201_ROADMAP_220331_tblLAB_RNA.csv')
dat <- read.csv(infile.rna,header=T)
dat$RNA_D <- as.Date(dat$RNA_D,format=c("%Y-%m-%d"))
dat <- data.table(dat)
# fit model/curve through each patient's VLs to get trajectory and interpolate
# set any values under/over threshold to one below/above the threshold
dat[RNA_V==-1 & !is.na(RNA_L), RNA_V:= RNA_L - 1]
dat[RNA_V==-1 & !is.na(RNA_UL), RNA_V:= RNA_UL + 1]
dat[RNA_V==-1, RNA_V:= 0]
dat[RNA_V==-1000 & RNA_L==-999, RNA_V:= -1]
# copy mean of each person's last 2 measurements to end of follow-up so we don't have missing values
dat2 <- dat %>%
  group_by(PATIENT) %>%
  slice_tail(n = 2)
dat2 <- data.table(dat2)
dat2 <- dat2[, list(RNA_V=mean(RNA_V,na.rm=T)),by='PATIENT']
dat2[, RNA_D:= max(dat$RNA_D)] # add date of last obs to impute until
dat <- merge(dat,dat2,by=c('PATIENT','RNA_D','RNA_V'),all=T)

# just keep patients in incidence cohort
dat <- subset(dat,PATIENT %in% dinf$SEQUENCE_ID)

# get infection date of recipient to predict VL of transmitter for

#dat <- merge(dat,subset(tmp,select=c('SEQUENCE_ID','EST_INF_DATE')),by.x=c('PATIENT'),by.y=c('SEQUENCE_ID'),all=T)
# add 1st jan for every year to predict
tmp <- data.table(expand.grid(PATIENT=unique(dat$PATIENT),
                  RNA_D=as.Date(c('2010-12-31','2011-12-31','2012-12-31','2013-12-31','2014-12-31','2015-12-31','2016-12-31',
                                  '2017-12-31','2018-12-31','2019-12-31','2020-12-31','2021-12-31'))))
dat <- merge(dat,tmp,by=c('PATIENT','RNA_D'),all=T)

# remove patients with less than 4 measurements
tmp <- unique(subset(dat,select=c('PATIENT','RNA_D','RNA_V')))
dn <- tmp[, list(N=length(RNA_V[!is.na(RNA_V)])),by=c('PATIENT')]
dat <- merge(dat,dn,by='PATIENT',all.x=T)
dat <- subset(dat, N>3)
dl <- subset(dat,PATIENT %notin% c('912511','904491','919461')) %>% # exclude  patients with problems fitting loess
  group_by(PATIENT) %>%
  arrange(PATIENT, RNA_D) %>%
  nest() %>%
  mutate(
    pred.response = purrr::map(data, function(x)stats::loess(RNA_V~as.numeric(RNA_D), span= 0.5, data = x,na.action="na.omit") %>%
                                 stats::predict(data.frame(RNA_D = as.numeric(x$RNA_D))))) %>%
  unnest(cols = c(data, pred.response))
dl <- data.table(dl)

# flag patients virally suppressed by putative infection date of recipient
dv <- data.table(unique(subset(dl,select=c('PATIENT','RNA_D','pred.response'))))
dv[pred.response<=200, supp:=1] # updated from 100 based on SHM report defn of suppression
dv[pred.response>200, supp:=0]
dv <- subset(dv,RNA_D %in% as.Date(c('2010-12-31','2011-12-31','2012-12-31','2013-12-31','2014-12-31','2015-12-31','2016-12-31',
                                       '2017-12-31','2018-12-31','2019-12-31','2020-12-31','2021-12-31')))
# add any pts back in who we could not predict VLs
tmp <- data.table(expand.grid(PATIENT=unique(dinf$SEQUENCE_ID),
                              RNA_D=as.Date(c('2010-12-31','2011-12-31','2012-12-31','2013-12-31','2014-12-31','2015-12-31','2016-12-31',
                                              '2017-12-31','2018-12-31','2019-12-31','2020-12-31','2021-12-31'))))
dv <- merge(dv,tmp,by=c('PATIENT','RNA_D'),all=T)

# merge in infection date
dv <- merge(dv,subset(dinf,select=c('SEQUENCE_ID','DIAGNOSIS_DATE','DIAGNOSIS_DATE_N','EST_INF_DATE','YEAR_OF_INF_EST','LOC_BIRTH_POS')),by.x='PATIENT',by.y='SEQUENCE_ID',all=T)

tmp <- dv[, list(#YEAR=c(2010,2010,2011,2011,2012,2012,2013,2013,2014,2014,2015,2015,2016,2016,2017,2017,
                 #       2018,2018,2019,2019,2020,2020,2021,2021),
                 N_inf_2010 = length(unique(PATIENT[EST_INF_DATE<='2010-12-31' & RNA_D>='2010-01-01' & RNA_D<='2010-12-31'])),
                 N_supp_2010 = length(unique(PATIENT[EST_INF_DATE<='2010-12-31' & RNA_D>='2010-01-01' & RNA_D<='2010-12-31' & supp==1])),
                 N_inf_2011 = length(unique(PATIENT[EST_INF_DATE<='2011-12-31' & RNA_D>='2011-01-01' & RNA_D<='2011-12-31'])),
                 N_supp_2011 = length(unique(PATIENT[EST_INF_DATE<='2011-12-31' & RNA_D>='2011-01-01' & RNA_D<='2011-12-31' & supp==1])),
                 N_inf_2012 = length(unique(PATIENT[EST_INF_DATE<='2012-12-31' & RNA_D>='2012-01-01' & RNA_D<='2012-12-31'])),
                 N_supp_2012 = length(unique(PATIENT[EST_INF_DATE<='2012-12-31' & RNA_D>='2012-01-01' & RNA_D<='2012-12-31' & supp==1])),
                 N_inf_2013 = length(unique(PATIENT[EST_INF_DATE<='2013-12-31' & RNA_D>='2013-01-01' & RNA_D<='2013-12-31'])),
                 N_supp_2013 = length(unique(PATIENT[EST_INF_DATE<='2013-12-31' & RNA_D>='2013-01-01' & RNA_D<='2013-12-31' & supp==1])),
                 N_inf_2014 = length(unique(PATIENT[EST_INF_DATE<='2014-12-31' & RNA_D>='2014-01-01' & RNA_D<='2014-12-31'])),
                 N_supp_2014 = length(unique(PATIENT[EST_INF_DATE<='2014-12-31' & RNA_D>='2014-01-01' & RNA_D<='2014-12-31' & supp==1])),
                 N_inf_2015 = length(unique(PATIENT[EST_INF_DATE<='2015-12-31' & RNA_D>='2015-01-01' & RNA_D<='2015-12-31'])),
                 N_supp_2015 = length(unique(PATIENT[EST_INF_DATE<='2015-12-31' & RNA_D>='2015-01-01' & RNA_D<='2015-12-31' & supp==1])),
                 N_inf_2016 = length(unique(PATIENT[EST_INF_DATE<='2016-12-31' & RNA_D>='2016-01-01' & RNA_D<='2016-12-31'])),
                 N_supp_2016 = length(unique(PATIENT[EST_INF_DATE<='2016-12-31' & RNA_D>='2016-01-01' & RNA_D<='2016-12-31' & supp==1])),
                 N_inf_2017 = length(unique(PATIENT[EST_INF_DATE<='2017-12-31' & RNA_D>='2017-01-01' & RNA_D<='2017-12-31'])),
                 N_supp_2017 = length(unique(PATIENT[EST_INF_DATE<='2017-12-31' & RNA_D>='2017-01-01' & RNA_D<='2017-12-31' & supp==1])),
                 N_inf_2018 = length(unique(PATIENT[EST_INF_DATE<='2018-12-31' & RNA_D>='2018-01-01' & RNA_D<='2018-12-31'])),
                 N_supp_2018 = length(unique(PATIENT[EST_INF_DATE<='2018-12-31' & RNA_D>='2018-01-01' & RNA_D<='2018-12-31' & supp==1])),
                 N_inf_2019 = length(unique(PATIENT[EST_INF_DATE<='2019-12-31' & RNA_D>='2019-01-01' & RNA_D<='2019-12-31'])),
                 N_supp_2019 = length(unique(PATIENT[EST_INF_DATE<='2019-12-31' & RNA_D>='2019-01-01' & RNA_D<='2019-12-31' & supp==1])),
                 N_inf_2020 = length(unique(PATIENT[EST_INF_DATE<='2020-12-31' & RNA_D>='2020-01-01' & RNA_D<='2020-12-31'])),
                 N_supp_2020 = length(unique(PATIENT[EST_INF_DATE<='2020-12-31' & RNA_D>='2020-01-01' & RNA_D<='2020-12-31' & supp==1])),
                 N_inf_2021 = length(unique(PATIENT[EST_INF_DATE<='2021-12-31' & RNA_D>='2021-01-01' & RNA_D<='2021-12-31'])),
                 N_supp_2021 = length(unique(PATIENT[EST_INF_DATE<='2021-12-31' & RNA_D>='2021-01-01' & RNA_D<='2021-12-31' & supp==1]))
),by='LOC_BIRTH_POS']
#tmp <- melt(tmp,id.vars='YEAR',measure.vars=patterns('^N_inf','^N_supp', cols=c('N_inf','N_supp')))
#tmp <- melt(tmp,id.vars='YEAR')
tmp2 <- melt(subset(tmp,select=c('LOC_BIRTH_POS',colnames(tmp)[grep('N_inf',colnames(tmp))])),id.vars='LOC_BIRTH_POS')
setnames(tmp2,c('variable','value'),c('YEAR','N_inf'))
tmp2[, YEAR:= gsub('N_inf_','',YEAR)]
#set(tmp2,NULL,'variable',NULL)
tmp3 <- melt(subset(tmp,select=c('LOC_BIRTH_POS',colnames(tmp)[grep('N_supp',colnames(tmp))])),id.vars='LOC_BIRTH_POS')
setnames(tmp3,c('variable','value'),c('YEAR','N_supp'))
tmp3[, YEAR:= gsub('N_supp_','',YEAR)]
#set(tmp3,NULL,'variable',NULL)
du <- merge(tmp2,tmp3,by=c('YEAR','LOC_BIRTH_POS'))

#tmp <- dcast(tmp,YEAR~c('N_inf','N_supp'),value.var = 'value')
#tmp[grep('N_inf',variable), N_inf:= value]
#tmp[grep('N_supp',variable), N_supp:= value]


# count unsuppressed
du[, unsupp:= N_inf - N_supp]

#tmp <- data.table(expand.grid(PATIENT=unique(dinf$SEQUENCE_ID),
#                              YEAR=as.Date(c('2010-12-31','2011-12-31','2012-12-31','2013-12-31','2014-12-31','2015-12-31','2016-12-31',
#                                              '2017-12-31','2018-12-31','2019-12-31','2020-12-31','2021-12-31'))))
#dv <- merge(dv,tmp,by=c('PATIENT'),all=T,allow.cartesian=T)
# count the number infected, and number suppressed by year
#du <- dv[, list(N_inf = length(unique(PATIENT[EST_INF_DATE<=YEAR])),
#                N_unsupp = sum(supp,na.rm=T[RNA_D==YEAR & EST_INF_DATE<=RNA_D])),by='YEAR']

#cat(paste0('Number of sources who were likely durably virally suppressed on infection date of recipient: ',nrow(subset(pairs,supp==1))))

#pairs <- subset(pairs,supp==0 | is.na(supp))
cat(" \n --------------------------------  calculate weighted average -------------------------------- \n")

# generate weights of total PLHIV across the 11 years
di <- dv[, list(N_inf=length(unique(PATIENT))),by=c('YEAR_OF_INF_EST','LOC_BIRTH_POS')]
dp <- du[, list(N_inf_tot=sum(N_inf)),by=c('LOC_BIRTH_POS')]
du <- merge(dp,du,by='LOC_BIRTH_POS')
du[, w:=N_inf/N_inf_tot]

wts <- tmp[LOC_BIRTH_POS=='Netherlands', list(YEAR=YEAR, w=N_total/sum(N_total))] # only do for NL because N_total same for all birth regions
tmp <- merge(tmp,wts,by='YEAR')
dp <- tmp[, list(pct=sum(prev*w)),by=c('LOC_BIRTH_POS')]
dp[, FROM_BPLACE:= factor(LOC_BIRTH_POS,
                          levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\n Caribbean','E. & C. Europe','MENA','Other'))]
saveRDS(dp,file=paste0(outfile.base,'-prevalence','.RDS'))

# calculate weights by region (% infections occurring in each year)

## calculate proportion of each region among the unsuppressed
# numerator = number unsuppressed (weighted av) across 2010-2021
# denominator = total unsuppressed (weighted av) across 2010-2021
