
## preamble ----
require(data.table)
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

outdir_undiag <- args$undiagnosed
job_tag_undiag <- args$job_tag_undiag

## read stanin
cat('\nReading Stan input data...')
infile.stanin <- list.files(args$outdir, pattern=paste0('_stanin.RData$'), recursive=TRUE)[1]
stopifnot(length(infile.stanin)>0)
stopifnot(length(infile.stanin)<=1)
tmp <- load(file.path(args$outdir, infile.stanin))
stopifnot(c('args','stan_data')%in%tmp)
args$analysis = 'analysis_220713'
args$indir = '~/Box\ Sync/Roadmap'
args$undiagnosed = '/Users/alexb/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Roadmap/sources/ethnicity_analysis/undiagnosed/undiagnosed_211102-cohort_2010_2015'

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

## calculate sampling prop for incident cases using diagnosed only ----

### load infection date data ----
dinf <- data.table(read.csv(file.path('data_Ams',args$analysis,'Infection_date_est_rec.csv')))
setnames(dinf,c("id",'estsctodiagMedian','estsctodiagLL','estsctodiagUL'),c("TO_SEQUENCE_ID",'SER_TO_DIAG','SER_TO_DIAG_LL','SER_TO_DIAG_UL'))
dinf <- unique(dinf)
dinf[,DIAGNOSIS_DATE:= as.Date(dinf[,hiv_pos_d],format="%Y-%m-%d")]
dinf[,DIAGNOSIS_DATE_N:= hivc.db.Date2numeric(dinf[,DIAGNOSIS_DATE])]
dinf[,EST_INF_DATE:= DIAGNOSIS_DATE_N-SER_TO_DIAG]
dinf[,EST_INF_DATE:= format(date_decimal(EST_INF_DATE), "%Y-%m-%d")]
dinf[,YEAR_OF_INF_EST := year(EST_INF_DATE)]

### merge in patient metadata ----
load(infile.meta)
dind <- data.table(dind)
dind[, SEQ:= PATIENT %in% ds$PATIENT]
dinf <- merge(dinf,subset(dind,select=c('PATIENT','CITY','TRANSM','LOC_BIRTH','ORIGIN','SEQ')),
              by.x='TO_SEQUENCE_ID',by.y='PATIENT',all.x=T)
dbas <- fread(infile.bas)
dbas[, DEATH_D := as.Date(DEATH_D,format="%Y-%m-%d")]
dinf <- merge(dinf,subset(dbas,select=c('PATIENT','DEATH_D')),by.x='TO_SEQUENCE_ID',by.y='PATIENT',all.x=T)

### add world regions ----
dinf[, LOC_BIRTH_POS:="Other"]
dinf[LOC_BIRTH %in% c("WEurope","NorthAm","Oceania"), LOC_BIRTH_POS:="W.Europe,\nN.America,Oceania"]
dinf[LOC_BIRTH %in% c("EEurope", "CEurope"), LOC_BIRTH_POS:="E. & C. Europe"]
dinf[LOC_BIRTH %in% c("LaAmCar"), LOC_BIRTH_POS:="S. America &\n Caribbean"]
dinf[LOC_BIRTH %in% c("DutchCarSuriname"), LOC_BIRTH_POS:="Suriname &\nDutch Caribbean"]
dinf[LOC_BIRTH %in% c("MENA"), LOC_BIRTH_POS:="MENA"]
dinf[ORIGIN=="NL", LOC_BIRTH_POS:="Netherlands"]

tmp <- subset(dinf,YEAR_OF_INF_EST>=2010)
sp <- tmp[, list(N=length(TO_SEQUENCE_ID),
                  N_seq = length(TO_SEQUENCE_ID[SEQ==T]),
                  pct=length(TO_SEQUENCE_ID)/nrow(tmp)),
           by='LOC_BIRTH_POS']

### merge in pct undiagnosed ----
ds <- readRDS(file=file.path(args$undiagnosed,paste0('p_undiagnosed_average_cohort_2010_2015_MSM.RDS')))
sp <- merge(sp,ds,by.x='LOC_BIRTH_POS',by.y='migrant_group')

### calculate total infected ----
sp[, N_inf:= N/(1-p0.5)]

### calculate sampling prob among incident cases ----
sp[, psi:= N_seq/N_inf]
sp[, N_inf_CL:=  N/(1-p0.025)]
sp[, N_inf_UL:=  N/(1-p0.975)]
sp[, psi_CL:=  N_seq/N_inf_CL]
sp[, psi_CU:=  N_seq/N_inf_UL]
saveRDS(sp,file=paste0(outfile.base,'-sampling_prob_cases','.RDS'))

# calculate number diagnosed from each birth region and number with a sequence

samples <- readRDS(file=file.path(outdir_undiag, paste0('samples_',job_tag_undiag,"_",args$trsm,'.rds')))
dmap <- readRDS(file=file.path(outdir, paste0("mapping_georeg_id.RDS")))

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

### calculate sp by year for cases ----
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

### calculate sp by year for sources ----
#mean_y <- merge(mean_y,dmap,by.x='mg',by.y='mgid')
mean_y <- readRDS(file=paste0(outfile.base,'-mean_undiagnosed_byyear_sources','.RDS'))
sp_s <- dinf[!is.na(YEAR_OF_INF_EST), list(N=length(TO_SEQUENCE_ID),
                    N_seq = length(TO_SEQUENCE_ID[SEQ==T])),
             by=c('LOC_BIRTH_POS','YEAR_OF_INF_EST')]
# add rows for years with missing values for each birth region? would need to make the N=0.0001 or similar..
sp_s <- merge(sp_s,mean_y,by.x=c('LOC_BIRTH_POS','YEAR_OF_INF_EST'),by.y=c('mwmb','year'),all.x=T)
sp_s[, N_inf:= N/(1-p0.5)]
sp_s[, psi:= N_seq/N_inf]
saveRDS(sp_s,file=paste0(outfile.base,'-sampling_prob_byyear_sources','.RDS')) # too many regions with missing years (since 1980) to do by year. Aggregate into 5/10yr bands?
sp_s <- sp_s[, list(N=sum(N),N_inf=sum(N_inf),N_seq=sum(N_seq)),by=c('LOC_BIRTH_POS')]
# calculate sampling prob among incident cases
sp_s[, psi:= N_seq/N_inf]
saveRDS(sp_s,file=paste0(outfile.base,'-sampling_prob_sources','.RDS'))


# calculate using birthplaces among incident cases in formulated pairs
#tmp2 <- unique(subset(do,select=c('TO_SEQUENCE_ID','TO_BPLACE')))
#tmp2 <- tmp2[, list(N_samp=length(unique(TO_SEQUENCE_ID)),
#                  pct_samp=length(unique(TO_SEQUENCE_ID))/nrow(tmp2)),
#           by='TO_BPLACE']

#sp <- merge(tmp,tmp2,by.x='LOC_BIRTH_POS',by.y='TO_BPLACE')
#sp[, psi:= N_samp/N]


## plot case-adjusted flows by birthplace ----
cat(" \n --------------------------------  plot adjusted flows by birthplace -------------------------------- \n")
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

g <- ggplot(subset(po,TO_BPLACE!='Overall')) + geom_bar(aes(x=TO_BPLACE,y=M,fill=FROM_BPLACE),stat='identity',position=position_dodge(width=0.9)) +
  geom_errorbar(aes(x=TO_BPLACE,ymin=CL, ymax=CU,fill=FROM_BPLACE),position=position_dodge(width=0.9), width=0.5, colour="black")	+
  scale_fill_npg() +
  labs(x='Birthplace of\nincident case',fill='Birthplace of likely\ntransmitter', y='Proportion of attributable infections\nto place of birth') +
  theme_bw(base_size=28) +
  theme(legend.pos='none',
        axis.title.x = element_blank(),
        #axis.text.x = element_text(angle=40, vjust = 0.5)) + #,
        axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9)) +#+ #,
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))

ggsave(file = paste0(outfile.base,'-adjusted_flowsINTO_samplingofcases_contributions.pdf'),
       g, w = 11, h = 8)
ggsave(file = paste0(outfile.base,'-adjusted_flowsINTO_samplingofcases_contributions.png'),
       g, w = 11, h = 8)


## adjust for sampling bias of sources ----
spy <- readRDS(file=paste0(outfile.base,'-sampling_prob_byyear_cases','.RDS'))
sp_s <- readRDS(file=paste0(outfile.base,'-sampling_prob_sources','.RDS')) # TODO: stratify estimates by 5 or 10yr intervals?

po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = 'tpair_prob' # should we be using tpair_prob now? unadjusted probs?
)
po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po <- data.table(po)
po[, PAIR_ID := as.integer(gsub(paste0('tpair_prob\\[([0-9]+)\\]'),'\\1',as.character(variable)))]
tmp <- subset(do, select = c('PAIR_ID','FROM_BPLACE','TO_BPLACE','FROM_SEQUENCE_ID','TO_SEQUENCE_ID'))
po <- merge(po, tmp, by = 'PAIR_ID')
setnames(po,'value','rho')

# expected number of missing trsm pairs m_j(z) from region z to case j
# predict counts from negbin using r = number of obs individuals per recipient, p = sampling prob
tmp2 <- subset(do,select=c('TO_SEQUENCE_ID','FROM_SEQUENCE_ID','FROM_BPLACE'))
tmp2 <- tmp2[, list(N=length(FROM_SEQUENCE_ID)),
             by=c('TO_SEQUENCE_ID','FROM_BPLACE')]

dmiss <- merge(subset(sp_s,select=c('LOC_BIRTH_POS','psi')),tmp2,by.x='LOC_BIRTH_POS',by.y='FROM_BPLACE')

# calculate median w(z) and multiply by the number of missing sources in cat z for recipient j (per MC sample)
dw <- po[, list(rho_med_src = median(rho)), by = c('draw','FROM_BPLACE')]
dw <- merge(dmiss,dw,by.x='LOC_BIRTH_POS',by.y='FROM_BPLACE',all.x=T,allow.cartesian = T) # adds rows per MC sample and per recipient missing obs for ALL world regions
# simulate missing sources per MC sample
#dw[, N_miss:= rnbinom(N,1,psi)]
#dw <- dw[, list(LOC_BIRTH_POS=LOC_BIRTH_POS, psi=psi,
#                 TO_SEQUENCE_ID= TO_SEQUENCE_ID,
#                 N=N,
#                 draw=draw,
#                 rho_med_src=rho_med_src,
#                 N_miss= rnbinom(n=nrow(dw),size=N,prob=psi))]
#setnames(dw,'LOC_BIRTH_POS','FROM_BPLACE')
# load imputed missing sources from poisson thinned model:
outdir <- "/Users/alexb/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Roadmap/sources/ethnicity_analysis/sampling_adjustments/case-specific-lambda"
dl <- readRDS(file=file.path(outdir, paste0('imputedsources_',job_tag,'.rds')))
setnames(dl,c('iter'),c('draw'))
dw <- merge(dw,dl,by=c('draw','TO_SEQUENCE_ID','LOC_BIRTH_POS'))
# sum the imputed rhos for the missing obs for each birthplace and case ID
mwz <- dw[, list(mwz = sum(N_miss * rho_med_src)),by=c('draw','TO_SEQUENCE_ID')]

#po <- merge(po, dw, by=c('draw','TO_SEQUENCE_ID','FROM_BPLACE')) # bug, dont merge sum of median missing to each pair

# calculate rho_kj = sum over observed rhos for incident case j (for all src cats/doesn't matter where from)
# sum over source cat z for unobserved sources for incident case j
tmp <- po[, list(rho_src = sum(rho)
                 #mwz = sum(N_miss * rho_med_src)
                 ),by=c('draw','TO_SEQUENCE_ID')]
po <- merge(po,tmp, by=c('draw','TO_SEQUENCE_ID'))
po <- merge(po,mwz, by=c('draw','TO_SEQUENCE_ID'))
po[rho_src + mwz==0, flag_zeroes:= 1]
dw <- merge(dw,tmp, by=c('draw','TO_SEQUENCE_ID'))
dw <- merge(dw,mwz, by=c('draw','TO_SEQUENCE_ID'))

# calculate new p_ij for obs pairs
# NB some values were NaN because rho, N_miss, rho_src and mwz are all 0. Set these to 0? or something small? 0.001?
dp_obs <- po[, list(p = rho/(rho_src + mwz)), by=c('draw','FROM_SEQUENCE_ID','TO_SEQUENCE_ID','FROM_BPLACE')]
dp_obs[is.nan(p), p:=0] # <<< TODO: check this

# calculate new p_ij for missing pairs (only for cases where we predict >0 missing sources)
dp_miss <- dw[N_miss>0, list(p = rho_med_src / (rho_src + mwz),
                     N_miss = N_miss), by=c('draw','FROM_BPLACE','TO_SEQUENCE_ID')]
#dp_miss[is.nan(p), p:=0] # <<< TODO: check this >>> UPDATE: no longer needed, no NaNs/infs
#dp_miss[is.infinite(p), p:=0] # <<< some are inf because denominator is 0 but numerator is non-zero (bc N_miss=0)

# calculate p_j(x) = sum p_ij for obs pairs over obs sources i from bplace z + sum p_j for missing pairs over # missing sources from bplace z
tmp <- dp_obs[, list(p_ijz = sum(p)),by=c('draw','TO_SEQUENCE_ID','FROM_BPLACE')]
tmp2 <- dp_miss[, list(p_jz = sum(p*N_miss)),by=c('draw','TO_SEQUENCE_ID','FROM_BPLACE')]

tmp <- merge(tmp,tmp2,by=c('draw','FROM_BPLACE','TO_SEQUENCE_ID'))
dp <- tmp[, list(pjx = p_ijz + p_jz), by=c('draw','FROM_BPLACE','TO_SEQUENCE_ID')] # sum to 1

# adjust for sampling of incident cases
dp <- merge(dp,unique(subset(do,select=c('TO_SEQUENCE_ID','TO_BPLACE','YEAR_OF_INF_EST'))),by='TO_SEQUENCE_ID')
dp <- merge(dp, subset(spy,select=c('LOC_BIRTH_POS','YEAR_OF_INF_EST','psi')),
            by.x=c('TO_BPLACE','YEAR_OF_INF_EST'),by.y=c('LOC_BIRTH_POS','YEAR_OF_INF_EST'))
# calculate flows from source cat Z
po <- dp[, list(value = sum(pjx/psi)), by = c('draw','FROM_BPLACE')] # adjust recipients for sampling too
tmp <- po[, list(total = sum(value)), by = c('draw')]
#tmp <- dp[, list(total = sum(pjx)), by = c('draw')]
po <- merge(po, tmp, by = 'draw')
po[, paf := value/total] # NOTE: we can actually just divide the sum of the pjx/psi by the number of unique recipients
saveRDS(po,file=paste0(outfile.base,'-flows_frombplace_adjusted_samplingbias_mcsamples','.RDS'))

po <- po[,
         list( q = quantile(paf, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('FROM_BPLACE')
]
po <- dcast.data.table(po, FROM_BPLACE~stat, value.var = 'q')
po[, TO_BPLACE:= 'Overall']
saveRDS(po,file=paste0(outfile.base,'-flows_frombplace_adjusted_samplingbias','.RDS'))


po[, FROM_BPLACE:= factor(FROM_BPLACE,
                          levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\n Caribbean','E. & C. Europe','MENA','Other'))]
pal <- pal_npg("nrc")(4)[c(1,3,4)]

g1 <- ggplot(subset(po,TO_BPLACE=='Overall')) +
  geom_bar(aes(x=FROM_BPLACE,y=M,fill=FROM_BPLACE),stat='identity',position=position_dodge(width=0.9)) +
  geom_errorbar(aes(x=FROM_BPLACE,ymin=CL, ymax=CU,fill=FROM_BPLACE),position=position_dodge(width=0.9),width=0.5, colour="black")	+
  scale_fill_npg() +
  labs(x='Birthplace of\nlikely source',fill='Birthplace of likely\ntransmitter', y='Proportion of attributable infections\nto place of birth') +
  theme_bw(base_size=28) +
  theme(legend.pos='none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9)) #+ #,

ggsave(file = paste0(outfile.base,'-adjusted_flows_frombplace_toplace_imputemissingdata.pdf'),
       g1, w = 11, h = 8)
ggsave(file = paste0(outfile.base,'-adjusted_flows_frombplace_toplace_imputemissingdata.png'),
       g1, w = 11, h = 8)


## stratified flows ----

# calculate flows from source cat Z to source cat Y
po <- dp[, list(value = sum(pjx/psi)), by = c('draw','FROM_BPLACE','TO_BPLACE')]
tmp <- po[, list(total = sum(value)), by = c('draw','TO_BPLACE')]
po <- merge(po, tmp, by = c('draw','TO_BPLACE'))
po[, paf := value/total]

po <- po[,
         list( q = quantile(paf, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('FROM_BPLACE','TO_BPLACE')
]
po <- dcast.data.table(po, TO_BPLACE + FROM_BPLACE ~stat, value.var = 'q')
saveRDS(po,file=paste0(outfile.base,'-stratified_flows_frombplace_adjusted_samplingbias','.RDS'))


po[, FROM_BPLACE:= factor(FROM_BPLACE,
                          levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\n Caribbean','E. & C. Europe','MENA','Other'))]
po[, TO_BPLACE:= factor(TO_BPLACE,
                          levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\n Caribbean','E. & C. Europe','MENA','Other'))]

g2 <- ggplot(po) + geom_bar(aes(x=TO_BPLACE,y=M,fill=FROM_BPLACE),stat='identity',position=position_dodge(width=0.9)) +
  geom_errorbar(aes(x=TO_BPLACE,ymin=CL, ymax=CU,fill=FROM_BPLACE),position=position_dodge(width=0.9), width=0.5, colour="black")	+
  scale_fill_npg() +
  labs(fill='Birthplace of\nlikely transmitter', y='Proportion of transmission flows',x='Birthplace of\nrecipient') +
  theme_bw() +
  theme(legend.pos='bottom') + #,
  #axis.text.x = element_text(angle=0, vjust = 0.5)) + #,
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))

ggsave(file = paste0(outfile.base,'-stratified_flows_frombplace_adjusted_samplingbias.pdf'),
       g2, w = 9, h = 6)
ggsave(file = paste0(outfile.base,'-stratified_flows_frombplace_adjusted_samplingbias.png'),
       g2, w = 9, h = 6)



g <- ggarrange(g1 + rremove("xlab") + theme_bw()+ theme(axis.text.x = element_text(angle=45, vjust = 0.9,hjust=0.9),legend.pos='bottom'),
               g2+ theme(axis.text.x = element_text(angle=45, vjust = 0.9,hjust=0.9))+ rremove("xlab") + rremove("ylab"),
               ncol=2,widths=c(0.35,0.65),align='h',common.legend=T,legend = "bottom")
#g_bplace <- annotate_figure(g,bottom = text_grob("Birthplace of recipient",size=28))

ggsave(file = paste0(outfile.base,'-flows_frombplace_adjusted_samplingbias_panel.pdf'),
       g, w = 9, h = 6)
ggsave(file = paste0(outfile.base,'-flows_frombplace_adjusted_samplingbias_panel.png'),
       g, w = 9, h = 6)

## plot case-adjusted flows by birthplace of case and source ----
cat(" \n --------------------------------  plot adjusted flows by birthplace -------------------------------- \n")
## % flows from each region PER region of birth of cases (i.e. sum to one per regino of birth of recipient)

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
po <- po[, list(value = sum(value/psi)), by = c('draw','FROM_BPLACE','TO_BPLACE')]
tmp <- po[, list(total = sum(value)), by = c('draw','TO_BPLACE')]
po <- merge(po, tmp, by = c('draw','TO_BPLACE'))
po[, paf := value/total]
po <- po[,
         list( q = quantile(paf, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('FROM_BPLACE','TO_BPLACE')
]
po <- dcast.data.table(po, FROM_BPLACE + TO_BPLACE ~stat, value.var = 'q')
po[, FROM_BPLACE:= factor(FROM_BPLACE,
                          levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\n Caribbean','E. & C. Europe','MENA','Other'))]
po[, TO_BPLACE:= factor(TO_BPLACE,
                        levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                 'S. America &\n Caribbean','E. & C. Europe','MENA','Other'))]
saveRDS(po,file=paste0(outfile.base,'-adjusted_flowsINTO_samplingofcases_bplacecase_bplacesrc','.RDS'))


## plot case-adjusted flows by birthplace of case and source ----
cat(" \n --------------------------------  plot adjusted flows by birthplace -------------------------------- \n")
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
po <- po[, list(value = sum(value/psi)), by = c('draw','FROM_BPLACE','TO_BPLACE')]
tmp <- po[, list(total = sum(value)), by = c('draw')]
po <- merge(po, tmp, by = 'draw')
po[, paf := value/total]
po <- po[,
         list( q = quantile(paf, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('FROM_BPLACE','TO_BPLACE')
]
po <- dcast.data.table(po, FROM_BPLACE + TO_BPLACE ~stat, value.var = 'q')
po[, FROM_BPLACE:= factor(FROM_BPLACE,
                          levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\n Caribbean','E. & C. Europe','MENA','Other'))]
po[, TO_BPLACE:= factor(TO_BPLACE,
                          levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\n Caribbean','E. & C. Europe','MENA','Other'))]
saveRDS(po,file=paste0(outfile.base,'-adjusted_flows_samplingofcases_bplacecase_bplacesrc','.RDS'))

