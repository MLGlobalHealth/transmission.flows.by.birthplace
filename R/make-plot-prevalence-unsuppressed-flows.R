
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
    job_tag = 'agegps_sensanalysis_210216_MSM-2010_2022',
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

cat(" \n --------------------------------  with arguments -------------------------------- \n")

infile.seq <-	file.path(args$indir, 'Data', 'data_220331/SHM_2201_ROADMAP_220331_tblLAB_seq.rda')
infile.bas <- file.path(args$indir, 'Data', 'data_220331','SHM_2201_ROADMAP_220331_tblBAS.csv')
infile.meta <- file.path(args$indir, args$analysis, 'misc', '220713_sequence_labels.rda')

## read stanin
cat('\nReading Stan input data...')
infile.stanin <- list.files(args$outdir, pattern=paste0('_stanin.RData$'), recursive=TRUE)[1]
stopifnot(length(infile.stanin)>0)
stopifnot(length(infile.stanin)<=1)
tmp <- load(file.path(args$outdir, infile.stanin))
stopifnot(c('args','stan_data')%in%tmp)
args$analysis = 'analysis_220713'
args$indir = '~/Box\ Sync/Roadmap'

tmp <- paste0(outfile.base,'-fitted_stan_model.rds')
cat("\n Read fitted dat ", tmp , "\n")
model_fit <- readRDS(file = tmp)

## load ethnicity data ----

load(infile.seq)

do <- merge(do,subset(dind,select=c('PATIENT','LOC_BIRTH')),by.x='FROM_SEQUENCE_ID',by.y='PATIENT',all.x=T)
setnames(do,'LOC_BIRTH','FROM_LOC_BIRTH')
do <- merge(do,subset(dind,select=c('PATIENT','LOC_BIRTH')),by.x='TO_SEQUENCE_ID',by.y='PATIENT',all.x=T)
setnames(do,'LOC_BIRTH','TO_LOC_BIRTH')

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

do[, FROM_MIGRANT:= 'Foreign-born']
do[FROM_ORIGIN=="NL", FROM_MIGRANT:= 'Dutch-born']
do[, TO_MIGRANT:= 'Foreign-born']
do[TO_ORIGIN=="NL", TO_MIGRANT:= 'Dutch-born']

do[, FROM_BPLACE:= factor(FROM_BPLACE,
                          levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\n Caribbean','E. & C. Europe','MENA','Other'))]
do[, TO_BPLACE:= factor(TO_BPLACE,
                        levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                 'S. America &\n Caribbean','E. & C. Europe','MENA','Other'))]

### calculate prevalence ----
## plot % birthplace among sources by % HIV positive in 2010-2021 (i.e. infection date <2021 or whole cohort)
## load infection date data
dinf <- data.table(read.csv(file.path('data_Ams',args$analysis,'Infection_date_est_rec.csv')))
setnames(dinf,c("id",'estsctodiagMedian','estsctodiagLL','estsctodiagUL'),c("TO_SEQUENCE_ID",'SER_TO_DIAG','SER_TO_DIAG_LL','SER_TO_DIAG_UL'))
dinf <- unique(dinf)
dinf[,DIAGNOSIS_DATE:= as.Date(dinf[,hiv_pos_d],format="%Y-%m-%d")]
dinf[,DIAGNOSIS_DATE_N:= hivc.db.Date2numeric(dinf[,DIAGNOSIS_DATE])]
dinf[,EST_INF_DATE:= DIAGNOSIS_DATE_N-SER_TO_DIAG]
dinf[,EST_INF_DATE:= format(date_decimal(EST_INF_DATE), "%Y-%m-%d")]
dinf[,YEAR_OF_INF_EST := year(EST_INF_DATE)]

## merge in patient metadata
load(infile.meta)
dinf <- merge(dinf,subset(dind,select=c('PATIENT','CITY','TRANSM','LOC_BIRTH','ORIGIN')),
              by.x='TO_SEQUENCE_ID',by.y='PATIENT',all.x=T)
dinf <- unique(subset(dinf,TRANSM=='MSM'))
dbas <- fread(infile.bas)
dbas[, DEATH_D := as.Date(DEATH_D,format="%Y-%m-%d")]
dinf <- merge(dinf,subset(dbas,select=c('PATIENT','DEATH_D')),by.x='TO_SEQUENCE_ID',by.y='PATIENT',all.x=T)

# calculate prevalence as proportion by birth region across years 2010-2021 (weighted average)
#tmp <- subset(dinf, (is.na(DEATH_D) | DEATH_D>='2010-01-01'))
dinf[, LOC_BIRTH_POS:="Other"]
dinf[LOC_BIRTH %in% c("WEurope","NorthAm","Oceania"), LOC_BIRTH_POS:="W.Europe,\nN.America,Oceania"]
dinf[LOC_BIRTH %in% c("EEurope", "CEurope"), LOC_BIRTH_POS:="E. & C. Europe"]
dinf[LOC_BIRTH %in% c("LaAmCar"), LOC_BIRTH_POS:="S. America &\n Caribbean"]
dinf[LOC_BIRTH %in% c("DutchCarSuriname"), LOC_BIRTH_POS:="Suriname &\nDutch Caribbean"]
dinf[LOC_BIRTH %in% c("MENA"), LOC_BIRTH_POS:="MENA"]
dinf[ORIGIN=="NL", LOC_BIRTH_POS:="Netherlands"]

# count the cases pre-2010
#pre <- dinf[,list(YEAR = 'pre-2010',
#                       N=length(TO_SEQUENCE_ID[YEAR_OF_INF_EST< 2010 & is.na(DEATH_D)])), # count total cases pre-2010, discounting anyone who died
#                 by=c('LOC_BIRTH_POS')]
#pre[, N_tot:= sum(N)]

# count the incident cases 2010-2021
tmp <- list()
tmp2 <- list()
for(i in 1980:2021){
  tmp[[i]] <- dinf[,list(YEAR = i,
                         N=length(TO_SEQUENCE_ID[YEAR_OF_INF_EST== i & is.na(DEATH_D)])),
                   by=c('LOC_BIRTH_POS')]
  tmp2[[i]] <- dinf[,list(YEAR = i,
                          N_total=length(TO_SEQUENCE_ID[YEAR_OF_INF_EST== i & is.na(DEATH_D)])),
  ]
}
tmp <- do.call(`rbind`,tmp)
tmp2 <- do.call(`rbind`,tmp2)
tmp <- merge(tmp,tmp2,by=c('YEAR'),all=T)

## adjust for undiagnosed ----
### load estimates of proportion undiagnosed

ds <- readRDS(file=file.path(args$undiagnosed,paste0('p_undiagnosed_byyear_cohort_2010_2015_MSM.RDS')))
ds <- dcast(ds,migrant_group+year~qlabel,value.var='p')
ds <- merge(tmp,ds,by.x=c('LOC_BIRTH_POS','YEAR'),by.y=c('migrant_group','year'),all=T)
ds[is.na(p0.5), p0.5:=0] # set pre-2010 prob(undiagnosed) to 0
## calculate total infected
ds[, N_inf:= round(N/(1-p0.5))]

## calculate prevalence 2010-2021
tmp <- list()
tmp2 <- list()
for(i in 1980:2021){
  tmp[[i]] <- ds[,list(YEAR = i,
                         N=sum(N_inf[YEAR<= i])),
                   by=c('LOC_BIRTH_POS')]
  tmp2[[i]] <- ds[,list(YEAR = i,
                          N_total=sum(N_inf[YEAR<= i])),
  ]
}
tmp <- do.call(`rbind`,tmp)
tmp2 <- do.call(`rbind`,tmp2)
tmp <- merge(tmp,tmp2,by=c('YEAR'),all=T)
tmp[, prev:= N/N_total]

# generate weights of total PLHIV across the 11 years
wts <- tmp[LOC_BIRTH_POS=='Netherlands' & YEAR %in% 2010:2021, list(YEAR=YEAR, w=N_total/sum(N_total))] # only do for NL because N_total same for all birth regions
tmp <- merge(tmp,wts,by='YEAR')
dp <- tmp[, list(pct=sum(prev*w)),by=c('LOC_BIRTH_POS')]
dp[, FROM_BPLACE:= factor(LOC_BIRTH_POS,
                          levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\n Caribbean','E. & C. Europe','MENA','Other'))]
saveRDS(dp,file=paste0(outfile.base,'-prevalence','.RDS'))

### plot prevalence ----
g_prev <- ggplot(subset(dp)) + geom_bar(aes(x=FROM_BPLACE,y=pct,fill=FROM_BPLACE),stat='identity',position=position_dodge(width=0.9)) +
  scale_fill_npg() +
  labs(x='Birthplace of\nincident case',fill='Birthplace of likely\ntransmitter', y='Prevalence by place of birth \n(weighted average 2010-2021)') +
  theme_bw(base_size=28) +
  theme(legend.pos='none',
        axis.title.x = element_blank(),
        #axis.text.x = element_text(angle=40, vjust = 0.5)) + #,
        axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9)) +#+ #,
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))
ggsave(file = paste0(outfile.base,'-rep_',replicate,'-contribution_to_prevalence.pdf'), g_prev, w = 23, h = 10)
ggsave(file = paste0(outfile.base,'-rep_',replicate,'-contribution_to_prevalence.png'), g_prev, w = 16, h = 10)

## plot flows by birthplace ----
cat(" \n --------------------------------  plot flows by birthplace -------------------------------- \n")
spy <- readRDS(file=paste0(outfile.base,'-sampling_prob_byyear_cases','.RDS'))

po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = 'tpair_prob_w' #tpair_prob_w
)
po <- data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po <- data.table(po)
po[, PAIR_ID := as.integer(gsub(paste0('tpair_prob_w\\[([0-9]+)\\]'),'\\1',as.character(variable)))]
tmp <- subset(do, select = c('PAIR_ID','FROM_BPLACE','TO_BPLACE','YEAR_OF_INF_EST'))
#setnames(tmp,'FROM_BPLACE','BPLACE')
po <- merge(po, tmp, by = 'PAIR_ID')
po <- merge(po, subset(spy,select=c('LOC_BIRTH_POS','YEAR_OF_INF_EST','psi')),
            by.x=c('TO_BPLACE','YEAR_OF_INF_EST'),by.y=c('LOC_BIRTH_POS','YEAR_OF_INF_EST'))
po <- po[, list(value = sum(value/psi)), by = c('draw','FROM_BPLACE')]
tmp <- po[, list(total = sum(value)), by = c('draw')]
po <- merge(po, tmp, by = 'draw')
po[, paf := value/total]
saveRDS(po,file=paste0(outfile.base,'-flows_frombplace_adjusted_samplingbias_mcsamples','.RDS'))

po <- po[,
         list( q = quantile(paf, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('FROM_BPLACE')
]
po <- dcast.data.table(po, FROM_BPLACE~stat, value.var = 'q')
setnames(po,'FROM_BPLACE','FROM_BPLACE')
po[, TO_BPLACE:= 'Overall']
saveRDS(po,file=paste0(outfile.base,'-rep_',replicate,'-adjusted_flows_samplingofcases','.RDS'))

po[, FROM_BPLACE:= factor(FROM_BPLACE,
                          levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\n Caribbean','E. & C. Europe','MENA','Other'))]

# for 3panel plot
g1 <- ggplot(subset(po,TO_BPLACE=='Overall')) + geom_bar(aes(x=FROM_BPLACE,y=M,fill=FROM_BPLACE),stat='identity',position=position_dodge(width=0.9)) +
  geom_errorbar(aes(x=FROM_BPLACE,ymin=CL, ymax=CU,fill=FROM_BPLACE),position=position_dodge(width=0.9),width=0.5, colour="black")	+
  scale_fill_npg() +
  labs(x='Birthplace of\nincident case',fill='Birthplace of likely\ntransmitter', y='Proportion of attributable infections\nto place of birth') +
  theme_bw(base_size=28) +
  theme(legend.pos='none',
        axis.title.x = element_blank(),
        #axis.text.x = element_text(angle=40, vjust = 0.5)) + #,
        axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9)) +#+ #,
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))
ggsave(file = paste0(outfile.base,'-rep_',replicate,'-contribution_to_flows.pdf'), g1, w = 23, h = 10)
ggsave(file = paste0(outfile.base,'-rep_',replicate,'-contribution_to_flows.png'), g1, w = 16, h = 10)

## plot flows over prevelance contributions ----
cat(" \n --------------------------------  plot flows over prevalence contributions -------------------------------- \n")

po <- readRDS(file=paste0(outfile.base,'-flows_frombplace_adjusted_samplingbias_mcsamples','.RDS'))

po <- merge(po, dp, by='FROM_BPLACE')
po[, contr_prev:= paf/pct]
po <- po[,
         list( q = quantile(contr_prev, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('FROM_BPLACE')
]
po <- dcast.data.table(po, FROM_BPLACE~stat, value.var = 'q')
po[, TO_BPLACE:= 'Overall']
saveRDS(po,file=paste0(outfile.base,'-rep_',replicate,'-flows_frombplace_hivpos','.RDS'))

po[, FROM_BPLACE:= factor(FROM_BPLACE,
                          levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\n Caribbean','E. & C. Europe','MENA','Other'))]
pal <- pal_npg("nrc")(4)[c(1,3,4)]

g2 <- ggplot(subset(po,TO_BPLACE=='Overall')) +
  geom_hline(yintercept=1,linetype=2) +
  geom_bar(aes(x=FROM_BPLACE,y=M,fill=FROM_BPLACE),stat='identity',position=position_dodge(width=0.9)) +
  geom_errorbar(aes(x=FROM_BPLACE,ymin=CL, ymax=CU,fill=FROM_BPLACE),position=position_dodge(width=0.9),width=0.5, colour="black")	+
  scale_fill_npg() +
  labs(x='Birthplace of\nincident case',fill='Birthplace of likely\ntransmitter', y='Contribution to transmission\nrelative to proportion of HIV positive') +
  theme_bw(base_size=28) +
  theme(legend.pos='none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9)) #+ #,
ggsave(file = paste0(outfile.base,'-rep_',replicate,'-contribution_to_flows_over_prev.pdf'), g2, w = 23, h = 10)
ggsave(file = paste0(outfile.base,'-rep_',replicate,'-contribution_to_flows_over_prev.png'), g2, w = 16, h = 10)

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
dinf <- unique(subset(dinf,YEAR_OF_INF_EST >= 2010 & YEAR_OF_INF_EST<2022))
dat <- subset(dat,PATIENT %in% dinf$TO_SEQUENCE_ID)

# get infection date of recipient to predict VL of transmitter for

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
tmp <- data.table(expand.grid(PATIENT=unique(dinf$TO_SEQUENCE_ID),
                              RNA_D=as.Date(c('2010-12-31','2011-12-31','2012-12-31','2013-12-31','2014-12-31','2015-12-31','2016-12-31',
                                              '2017-12-31','2018-12-31','2019-12-31','2020-12-31','2021-12-31'))))
dv <- merge(dv,tmp,by=c('PATIENT','RNA_D'),all=T)

# merge in infection date
dv <- merge(dv,subset(dinf,select=c('TO_SEQUENCE_ID','DIAGNOSIS_DATE','DIAGNOSIS_DATE_N','EST_INF_DATE','YEAR_OF_INF_EST','LOC_BIRTH_POS')),
            by.x='PATIENT',by.y='TO_SEQUENCE_ID',all=T)

tmp <- dv[, list(N_inf_2010 = length(unique(PATIENT[EST_INF_DATE<='2010-12-31' & RNA_D>='2010-01-01' & RNA_D<='2010-12-31'])),
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
tmp2 <- melt(subset(tmp,select=c('LOC_BIRTH_POS',colnames(tmp)[grep('N_inf',colnames(tmp))])),id.vars='LOC_BIRTH_POS')
setnames(tmp2,c('variable','value'),c('YEAR','N_inf'))
tmp2[, YEAR:= gsub('N_inf_','',YEAR)]
tmp3 <- melt(subset(tmp,select=c('LOC_BIRTH_POS',colnames(tmp)[grep('N_supp',colnames(tmp))])),id.vars='LOC_BIRTH_POS')
setnames(tmp3,c('variable','value'),c('YEAR','N_supp'))
tmp3[, YEAR:= gsub('N_supp_','',YEAR)]
du <- merge(tmp2,tmp3,by=c('YEAR','LOC_BIRTH_POS'))

# count unsuppressed
du[, unsupp:= N_inf - N_supp]
set(du,NULL,'N_inf',NULL)

tmp <- du[, list(unsupp_tot=sum(unsupp)),by=c('YEAR')]
du <- merge(du,tmp,by='YEAR')

## calculate weighted average ----
cat(" \n --------------------------------  calculate weighted average -------------------------------- \n")

# generate weights of total PLHIV across the 11 years
#di <- dv[, list(N_inf=length(unique(PATIENT))),by=c('YEAR_OF_INF_EST','LOC_BIRTH_POS')]
#dp <- sp[, list(N_inf_tot=sum(N_inf)),by=c('LOC_BIRTH_POS')]
#sp <- merge(sp,dp,by='LOC_BIRTH_POS')
sp[, w:=N_inf/N_inf_tot]

du[, YEAR:= as.integer(YEAR)]
du <- merge(du,sp,by.x=c('LOC_BIRTH_POS','YEAR'),by.y=c('LOC_BIRTH_POS','YEAR_OF_INF_EST'),all=T)
du[is.na(w), w:=0]
du[is.na(N_inf), N_inf:=0]
du[is.na(N_inf_tot), N_inf_tot:=0]

# calculate proportion of birth regions among the unsuppressed as a weighted sum of the years
dp <- du[, list(p_among_unsupp=sum((unsupp/unsupp_tot)*w)),by=c('LOC_BIRTH_POS')]

dp[, FROM_BPLACE:= factor(LOC_BIRTH_POS,
                          levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\n Caribbean','E. & C. Europe','MENA','Other'))]
saveRDS(dp,file=paste0(outfile.base,'-bplace_unsuppressed_among_PLHIV','.RDS'))

## calculate flows/unsuppressed
cat(" \n --------------------------------  calculate flows/unsuppressed -------------------------------- \n")
spy <- readRDS(file=paste0(outfile.base,'-sampling_prob_byyear_cases','.RDS'))

po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = 'tpair_prob_w' #tpair_prob_w
)
po <- data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po <- data.table(po)
po[, PAIR_ID := as.integer(gsub(paste0('tpair_prob_w\\[([0-9]+)\\]'),'\\1',as.character(variable)))]
tmp <- subset(do, select = c('PAIR_ID','FROM_BPLACE','TO_BPLACE','YEAR_OF_INF_EST'))
#setnames(tmp,'FROM_BPLACE','BPLACE')
po <- merge(po, tmp, by = 'PAIR_ID')
po <- merge(po, subset(spy,select=c('LOC_BIRTH_POS','YEAR_OF_INF_EST','psi')),
            by.x=c('TO_BPLACE','YEAR_OF_INF_EST'),by.y=c('LOC_BIRTH_POS','YEAR_OF_INF_EST'))
po <- po[, list(value = sum(value/psi)), by = c('draw','FROM_BPLACE')]
tmp <- po[, list(total = sum(value)), by = c('draw')]
po <- merge(po, tmp, by = 'draw')
po[, paf := value/total]

saveRDS(po,file=paste0(outfile.base,'-adjusted_flows_samplingofcases_MCsamples','.RDS'))

po <- merge(po,dp,by=c('FROM_BPLACE'))
po[, contr_unsupp:= paf/p_among_unsupp]
po <- po[,
         list( q = quantile(contr_unsupp, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('FROM_BPLACE')
]
po <- dcast.data.table(po, FROM_BPLACE~stat, value.var = 'q')
setnames(po,'FROM_BPLACE','FROM_BPLACE')
po[, TO_BPLACE:= 'Overall']

po[, FROM_BPLACE:= factor(FROM_BPLACE,
                          levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\n Caribbean','E. & C. Europe','MENA','Other'))]

saveRDS(po,file=paste0(outfile.base,'-adjusted_flows_rel_unsupp_samplingofcases','.RDS'))


## plot unsuppressed ----
g_uns <- ggplot(subset(dp)) + geom_bar(aes(x=FROM_BPLACE,y=p_among_unsupp,fill=FROM_BPLACE),stat='identity',position=position_dodge(width=0.9)) +
  scale_fill_npg() +
  labs(x='Birthplace of\nincident case',fill='Birthplace of likely\ntransmitter', y='Composition of unsuppressed among PLHIV \nby place of birth \n(weighted average 2010-2021)') +
  theme_bw(base_size=28) +
  theme(legend.pos='none',
        axis.title.x = element_blank(),
        #axis.text.x = element_text(angle=40, vjust = 0.5)) + #,
        axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9)) +#+ #,
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))
ggsave(file = paste0(outfile.base,'-rep_',replicate,'-contribution_to_unsuppressed.pdf'), g_uns, w = 23, h = 10)
ggsave(file = paste0(outfile.base,'-rep_',replicate,'-contribution_to_unsuppressed.png'), g_uns, w = 16, h = 10)

## plot ratio unsuppressed ----

pal <- pal_npg("nrc")(4)[c(1,3,4)]

g3 <- ggplot(subset(po,TO_BPLACE=='Overall')) +
  geom_hline(yintercept=1,linetype=2) +
  geom_bar(aes(x=FROM_BPLACE,y=M,fill=FROM_BPLACE),stat='identity',position=position_dodge(width=0.9)) +
  geom_errorbar(aes(x=FROM_BPLACE,ymin=CL, ymax=CU,fill=FROM_BPLACE),position=position_dodge(width=0.9),width=0.5, colour="black")	+
  scale_fill_npg() +
  labs(x='Birthplace of\nincident case',fill='Birthplace of likely\ntransmitter', y='Contribution to transmissions\nrelative to contribution towards PLHIV\nwith an unsuppressed virus') +
  theme_bw(base_size=28) +
  theme(legend.pos='none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9)) #+ #,
ggsave(file = paste0(outfile.base,'-rep_',replicate,'-flows_frombplace_contributions_rel_unsuppressed_contributions.pdf'), g_ratio, w = 23, h = 10)
ggsave(file = paste0(outfile.base,'-rep_',replicate,'-flows_frombplace_contributions_rel_unsuppressed_contributions.png'), g_ratio, w = 16, h = 10)

## Make panel plot ----

g <- ggarrange(g_prev+theme_bw(base_size=11) + labs(y='\nContribution to prevalence\n2010-2021') + theme(legend.pos='none',
                                                         axis.title.x = element_blank(),axis.text.x = element_text(angle=50, vjust = 0.95,hjust = 0.9)),
                      g_uns+theme_bw(base_size=11) + labs(y='\nContribution to PLHIV\nwith an unsuppressed virus') + theme(legend.pos='none',
                                                           axis.title.x = element_blank(),axis.text.x = element_text(angle=50, vjust = 0.95,hjust = 0.9)),
                      g1+theme_bw(base_size=11) + labs(y='\nContribution to transmissions') + theme(legend.pos='none',
                                                         axis.title.x = element_blank(),axis.text.x = element_text(angle=50, vjust = 0.95,hjust = 0.9)),
                      g2+theme_bw(base_size=11) + labs(y='\nContribution to transmissions\n relative to contribution to prevalence\n2010-2021') + theme(legend.pos='none',
                                                        axis.title.x = element_blank(),axis.text.x = element_text(angle=50, vjust = 0.95,hjust = 0.9)),
                      g3+theme_bw(base_size=11) + labs(y='\nContribution to transmissions\nrelative to contribution towards PLHIV\nwith an unsuppressed virus') + theme(legend.pos='none',
                                                        axis.title.x = element_blank(),axis.text.x = element_text(angle=50, vjust = 0.95,hjust = 0.9)),
                      ncol=3,nrow=2,align='v',labels='AUTO',font.label=list(size=14))

g_bplace <- annotate_figure(g,bottom = text_grob("Birthplace of incident case",size=11))
ggsave(file = paste0(outfile.base,'-prev_uns_flows_frombplace_contributions_panel.pdf'),
       g_bplace, w = 11, h = 7)
ggsave(file = paste0(outfile.base,'-prev_uns_flows_frombplace_contributions_panel.png'),
       g_bplace, w = 11, h = 7)

