
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
require(lubridate)

if (1)
{
  args <- list(
    source_dir = '~/Documents/GitHub/transmission.flows.by.birthplace',
    indir = '~/Box\ Sync/Roadmap/source_attribution',
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

## read stanin
cat('\nReading Stan input data...')
infile.stanin <- list.files(args$outdir, pattern=paste0('_stanin.RData$'), recursive=TRUE)[1]
stopifnot(length(infile.stanin)>0)
stopifnot(length(infile.stanin)<=1)
tmp <- load(file.path(args$outdir, infile.stanin))
stopifnot(c('args','stan_data')%in%tmp)

tmp <- paste0(outfile.base,'-fitted_stan_model.rds')
cat("\n Read fitted dat ", tmp , "\n")
model_fit <- readRDS(file = tmp)

## load ethnicity data ----

infile.seq <-	file.path(args$indir, 'Data', 'data_220331/SHM_2201_ROADMAP_220331_tblLAB_seq.rda')
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

# define stage in cascade of source at time of infection ----

#do[FROM_DIAGNOSIS_DATE>TO_EST_INFECTION_DATE, stage:= 'undiagnosed']
do[FROM_EST_INFECTION_DATE + years(1) <= TO_EST_INFECTION_DATE & FROM_DIAGNOSIS_DATE>TO_EST_INFECTION_DATE, STAGE:= 'undiag_first_year']
do[FROM_EST_INFECTION_DATE + years(1) > TO_EST_INFECTION_DATE & FROM_DIAGNOSIS_DATE>TO_EST_INFECTION_DATE, STAGE:= 'undiag_more_1_year']
do[FROM_DIAGNOSIS_DATE<=TO_EST_INFECTION_DATE & (supp==0 | is.na(supp)), STAGE:= 'diagnosed_unsuppressed']# <= because 3 were NA - because they were diagnosed on the estimated infection date of their recipient
do[FROM_DIAGNOSIS_DATE<=TO_EST_INFECTION_DATE & supp==1, STAGE:= 'suppressed']# <= because 3 were NA - because they were diagnosed on the estimated infection date of their recipient


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

tmp2 <- list()
for(i in 2010:2021){
  tmp2[[i]] <- dinf[,list(YEAR = i,
                          N_total=length(TO_SEQUENCE_ID[YEAR_OF_INF_EST<= i & (is.na(DEATH_D) | year(DEATH_D)>=YEAR_OF_INF_EST)])),
  ]
}
tmp <- do.call(`rbind`,tmp)
tmp2 <- do.call(`rbind`,tmp2)

### calculate % unsuppressed among population as weighted average ----
# load VL data

infile.rna <-	file.path(args$indir, 'Data', 'data_220331/SHM_2201_ROADMAP_220331_tblLAB_RNA.csv')
dat <- fread(file=infile.rna,header=T)
dat$RNA_D <- as.Date(dat$RNA_D,format=c("%Y-%m-%d"))
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

# remove patients with less than 2 measurements
tmp <- unique(subset(dat,select=c('PATIENT','RNA_D','RNA_V')))
dn <- tmp[, list(N=length(RNA_V[!is.na(RNA_V) & RNA_V>0])),by=c('PATIENT')]
dat <- merge(dat,dn,by='PATIENT',all.x=T)
dat <- subset(dat, N>3)
dl <- dat %>%
  tidyr::nest(-PATIENT) %>%
  dplyr::mutate(
    # Perform loess calculation on each PATIENT - fit to log(VL) and exp predictions to avoid zeroes
    m = purrr::map(data, loess,
                   formula = log(RNA_V + 0.1) ~ as.numeric(RNA_D), span = .6,na.action="na.omit"),
    # Retrieve the fitted values from each model
    fitted = purrr::map(m, `[[`, "fitted")
  )

names(dl$m) <-  dl$PATIENT
new <- purrr::map(dl$m, ~ predict(.x, newdata = as.Date(c('2010-01-01','2011-01-01','2012-01-01','2013-01-01',
                                                          '2014-01-01','2015-01-01','2016-01-01','2017-01-01',
                                                          '2018-01-01','2019-01-01','2020-01-01','2021-01-01'))))
new2 <- do.call(`rbind`,new)
rownames(new2) <- dl$PATIENT
new2 <- data.table(reshape2::melt(new2))
setnames(new2,c('Var1','Var2','value'),c('PATIENT','ID','pred_VL'))
new2[, pred_VL:= exp(pred_VL)]
new2 <- new2[order(PATIENT,ID),]
tmp <- data.table(ID=seq(1,12,1),
                  YEAR=seq(2010,2021,1))
new2 <- merge(new2,tmp,by='ID')
new2 <- new2[order(PATIENT,YEAR),]
new2 <- merge(new2,tmp2,by='YEAR')

# calculate proportion unsuppressed at start of each year by world region
# NB there are some missing individuals without infection dates but with VLs - why?
new2 <- merge(new2,subset(dinf,select=c('TO_SEQUENCE_ID','LOC_BIRTH_POS','DIAGNOSIS_DATE_N','EST_INF_DATE')),by.x='PATIENT',by.y='TO_SEQUENCE_ID',all.x=T)
#new2 <- subset(new2, !is.na(LOC_BIRTH_POS))

tmp <- new2[,list(N_unsupp=length(PATIENT[pred_VL>200 & !is.na(pred_VL)])),
            by=c('LOC_BIRTH_POS','YEAR')]
tmp2 <- tmp[, list(N_total=sum(N)),by='YEAR']
tmp <- merge(tmp, tmp2, by='YEAR')
tmp[, pct_uns:= N/N_total]

# generate weights of total PLHIV across the 11 years
tmp <- merge(tmp,wts,by='YEAR')
dp <- tmp[, list(pct=sum(pct_uns*w)),by=c('LOC_BIRTH_POS')]

dp[, FROM_BPLACE:= factor(LOC_BIRTH_POS,
                          levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\n Caribbean','E. & C. Europe','MENA','Other'))]
saveRDS(dp,file=paste0(outfile.base,'-unsuppresssed','.RDS'))

## classify sources by cascade stage ----
do[FROM_EST_INFECTION_DATE + years(1) <= TO_EST_INFECTION_DATE & FROM_DIAGNOSIS_DATE>TO_EST_INFECTION_DATE, STAGE:= 'undiag_first_year']
do[FROM_EST_INFECTION_DATE + years(1) > TO_EST_INFECTION_DATE & FROM_DIAGNOSIS_DATE>TO_EST_INFECTION_DATE, STAGE:= 'undiag_more_1_year']
do[FROM_DIAGNOSIS_DATE<=TO_EST_INFECTION_DATE & (supp==0 | is.na(supp)), STAGE:= 'diagnosed_unsuppressed']# <= because 3 were NA - because they were diagnosed on the estimated infection date of their recipient
do[FROM_DIAGNOSIS_DATE<=TO_EST_INFECTION_DATE & supp==1, STAGE:= 'suppressed']# <= because 3 were NA - because they were diagnosed on the estimated infection date of their recipient

### calculate sampling fraction ----

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

### calculate sampling prob ----
sp[, psi:= N_seq/N]
sp[, psi_und:= N_seq/N_inf]

## calculate sampling frac overall

spall <- sp[, list(N=sum(N),
                 N_seq = sum(N_seq),
                 N_inf = sum(N_inf))]

### calculate sampling prob ----
spall[, psi:= N_seq/N]
spall[, psi_und:= N_seq/N_inf]


### sampling prob undiagnosed ----
samples <- readRDS(file=file.path(outdir, paste0('samples_',job_tag,"_",args$trsm,'.rds')))
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

dat <- tidyr::crossing(year=seq(2010,2021,1),month=seq(1,12,1))
#dat <- tidyr::crossing(year=seq(2010,2021,1))
ds <- merge(dat,ds,all=T)
ds <- data.table(ds)
ds[, time:=(2022+(1/12))-(year + (month/12))] # +(1/12) to add 1 month so we count until end Dec 2018/start of Jan 2019. but I think it makes sense someone infected Dec 2018 has 0 prob of being diagnosed by end of 2018 due to time taken to detect virus
#ds[, time:=(2022-year)] # +(1/12) to add 1 month so we count until end Dec 2018/start of Jan 2019. but I think it makes sense someone infected Dec 2018 has 0 prob of being diagnosed by end of 2018 due to time taken to detect virus
ds[, time_1yr:= 1] # +(1/12) to add 1 month so we count until end Dec 2018/start of Jan 2019. but I think it makes sense someone infected Dec 2018 has 0 prob of being diagnosed by end of 2018 due to time taken to detect virus
ds[year==2021, time_1yr:= 1+(1/12) - (month/12)] # +(1/12) to add 1 month so we count until end Dec 2018/start of Jan 2019. but I think it makes sense someone infected Dec 2018 has 0 prob of being diagnosed by end of 2018 due to time taken to detect virus
ds[, p:=1 - pweibull(time,shape=shape,scale=scale)]
ds[, p_1yr:=1 - pweibull(time_1yr,shape=shape,scale=scale)]

mean_y <- ds[, list(p=quantile(p,prob=c(0.025,0.5,0.975)),
                    p_1yr=quantile(p_1yr,prob=c(0.025,0.5,0.975)),
                    qlabel=c('p0.025','p0.5','p0.975')),
             by=c('trsm','mg','year')] # summarise quantiles for each year
mean_y <- merge(mean_y,dmap,by.x='mg',by.y='mgid')
mean_y <- dcast(mean_y,trsm+mwmb+year~qlabel,value.var=c("p","p_1yr"))

tmp <- subset(dinf,YEAR_OF_INF_EST>=2010)
spu <- tmp[, list(N=length(TO_SEQUENCE_ID),
                  N_seq = length(TO_SEQUENCE_ID[SEQ==T]),
                  pct=length(TO_SEQUENCE_ID)/nrow(tmp)),
           by=c('LOC_BIRTH_POS','YEAR_OF_INF_EST')]
spu <- merge(spu,mean_y,by.x=c('LOC_BIRTH_POS','YEAR_OF_INF_EST'),by.y=c('mwmb','year'))
spu[, N_inf:= N/(1-p_p0.5)]
spu[, N_inf_1yr:= N/(1-p_1yr_p0.5)] ## need to revise this
spu[, psi_u:= N_seq/N_inf]

## plot posterior prob of being a pair by birthplace and stage in cascade ----
cat(" \n --------------------------------  plot flows by birthplace -------------------------------- \n")

po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = 'tpair_prob_w'
)
po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po <- data.table(po)
po[, PAIR_ID := as.integer(gsub(paste0('tpair_prob_w\\[([0-9]+)\\]'),'\\1',as.character(variable)))]
po <- po[,
         list( q = quantile(value, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('PAIR_ID')
]
po <- dcast.data.table(po, PAIR_ID~stat, value.var = 'q')
tmp <- subset(do, select = c('PAIR_ID','TO_BPLACE','STAGE','TIME_ELAPSED'))
po <- merge(po, tmp, by = 'PAIR_ID')

saveRDS(po,file=paste0(outfile.base,'-tpair_tobplace_stage','.RDS'))

po[, STAGE:= factor(STAGE,
                    levels=c('undiag_first_year','undiag_more_1_year','diagnosed_unsuppressed'),
                    labels=c('Undiagnosed (first year of infection)',
                             'Undiagnosed (after first year of infection)',
                             'Diagnosed but not virally suppressed'))]

p <- ggplot(po, aes(x = TIME_ELAPSED, color = TO_BPLACE)) +
  geom_point(aes(y = M)) +
  geom_errorbar(aes(ymin = IL, ymax = IU)) +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  ggsci::scale_color_npg() +
  labs(x = 'time elapsed (years)', y = 'transmission pair probability', colour = 'source category of transmitter') +
  facet_grid(~STAGE, scales = "free_x") +
  theme_bw() +
  theme( legend.position = 'bottom' )
ggsave(file = paste0(outfile.base,'-pi_ij_timeelapsed_tobplace_stage.png'), p, w = 12, h = 7)
ggsave(file = paste0(outfile.base,'-pi_ij_timeelapsed_tobplace_stage.pdf'), p, w = 12, h = 7)


po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = 'tpair_prob_w'
)
po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po[, PAIR_ID := as.integer(gsub(paste0('tpair_prob_w\\[([0-9]+)\\]'),'\\1',as.character(variable)))]
po <- po[,
         list( q = quantile(value, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('PAIR_ID')
]
po <- dcast.data.table(po, PAIR_ID~stat, value.var = 'q')
tmp <- subset(do, select = c('PAIR_ID','FROM_BPLACE','STAGE','TIME_ELAPSED'))
po <- merge(po, tmp, by = 'PAIR_ID')

saveRDS(po,file=paste0(outfile.base,'-tpair_frombplace_stage','.RDS'))

po[, STAGE:= factor(STAGE,
                    levels=c('undiag_first_year','undiag_more_1_year','diagnosed_unsuppressed'),
                    labels=c('Undiagnosed (first year of infection)',
                             'Undiagnosed (after first year of infection)',
                             'Diagnosed but not virally suppressed'))]

p <- ggplot(po, aes(x = TIME_ELAPSED, color = FROM_BPLACE)) +
  geom_point(aes(y = M)) +
  geom_errorbar(aes(ymin = IL, ymax = IU)) +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  ggsci::scale_color_npg() +
  labs(x = 'time elapsed (years)', y = 'transmission pair probability', colour = 'source category of transmitter') +
  facet_grid(~STAGE, scales = "free_x") +
  theme_bw() +
  theme( legend.position = 'bottom' )
ggsave(file = paste0(outfile.base,'-pi_ij_timeelapsed_frombplace_stage.png'), p, w = 12, h = 7)
ggsave(file = paste0(outfile.base,'-pi_ij_timeelapsed_frombplace_stage.pdf'), p, w = 12, h = 7)


## plot flows by birthplace and stage in cascade ----
cat(" \n --------------------------------  plot flows by birthplace -------------------------------- \n")

po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = 'tpair_prob_w'
)
po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po <- data.table(po)
po[, PAIR_ID := as.integer(gsub(paste0('tpair_prob_w\\[([0-9]+)\\]'),'\\1',as.character(variable)))]
tmp <- subset(do, select = c('PAIR_ID','TO_BPLACE','STAGE'))
po <- merge(po, tmp, by = 'PAIR_ID')
po <- merge(po, subset(sp,select=c('LOC_BIRTH_POS','psi')),by.x='TO_BPLACE',by.y='LOC_BIRTH_POS')
po <- po[, list(value = sum(value/psi)), by = c('draw','TO_BPLACE','STAGE')] # adjust recipients for sampling
#po <- po[, list(value = sum(value)), by = c('draw','TO_BPLACE','STAGE')]
tmp <- po[, list(total = sum(value)), by = c('draw','TO_BPLACE')]
po <- merge(po, tmp, by = c('draw','TO_BPLACE'))
po[, paf := value/total]

po <- po[,
         list( q = quantile(paf, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('TO_BPLACE','STAGE')
]
po <- dcast.data.table(po, TO_BPLACE+STAGE~stat, value.var = 'q')
saveRDS(po,file=paste0(outfile.base,'-flows_tobplace_stage','.RDS'))

po[, STAGE:= factor(STAGE,
                    levels=c('undiag_first_year','undiag_more_1_year','diagnosed_unsuppressed'),
                    labels=c('Undiagnosed (first year of infection)',
                             'Undiagnosed (after first year of infection)',
                             'Diagnosed but not virally suppressed'))]
g <- ggplot(subset(po)) + geom_bar(aes(fill=STAGE,y=M,x=TO_BPLACE),stat='identity',position=position_dodge(width=0.9)) +
  geom_errorbar(aes(fill=STAGE,ymin=CL, ymax=CU,x=TO_BPLACE),position=position_dodge(width=0.9), width=0.5, colour="black")	+
  scale_fill_npg('nrc') +
  #scale_fill_manual(name="Stage in cascade of\nlikely transmitter at\nputative infection time",values = c('Dutch-born'=pal[2],'Foreign-born'=pal[3])) +
  labs(x='Birthplace of recipient', y='Proportion of attributable\ninfections to place of birth', fill = 'Stage in cascade of likely transmitter at putative infection time') +
  theme_bw(base_size=28) +
  theme(legend.pos='bottom',
        #axis.text.y=element_blank(),
        strip.background=element_blank(),
        strip.text=element_blank(),
        axis.text.x = element_text(angle=45,vjust=0.7)) + #,
  guides(fill = guide_legend(nrow = 3,title.position='top'),by.col=T) +
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))
ggsave(file = paste0(outfile.base,'-flows_tobplace_stage.pdf'), g , w = 16, h = 12)
ggsave(file = paste0(outfile.base,'-flows_tobplace_stage.png'), g , w = 16, h = 12)


## plot flows by birthplace of source and stage in cascade ----
cat(" \n --------------------------------  plot flows by birthplace of source -------------------------------- \n")

po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = 'tpair_prob_w'
)
po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po <- data.table(po)
po[, PAIR_ID := as.integer(gsub(paste0('tpair_prob_w\\[([0-9]+)\\]'),'\\1',as.character(variable)))]
tmp <- subset(do, select = c('PAIR_ID','FROM_BPLACE','STAGE'))
po <- merge(po, tmp, by = 'PAIR_ID')
po <- merge(po, subset(sp,select=c('STAGE','psi')),by.x='STAGE',by.y='STAGE')
po <- po[, list(value = sum(value/psi)), by = c('draw','FROM_BPLACE','STAGE')] # adjust recipients for sampling
#po <- po[, list(value = sum(value)), by = c('draw','FROM_BPLACE','STAGE')]
tmp <- po[, list(total = sum(value)), by = c('draw','FROM_BPLACE')]
po <- merge(po, tmp, by = c('draw','FROM_BPLACE'))
po[, paf := value/total]

po <- po[,
         list( q = quantile(paf, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('FROM_BPLACE','STAGE')
]
po <- dcast.data.table(po, FROM_BPLACE+STAGE~stat, value.var = 'q')
saveRDS(po,file=paste0(outfile.base,'-flows_frombplace_stage','.RDS'))

po[, STAGE:= factor(STAGE,
                    levels=c('undiag_first_year','undiag_more_1_year','diagnosed_unsuppressed'),
                    labels=c('Undiagnosed (first year of infection)',
                             'Undiagnosed (after first year of infection)',
                             'Diagnosed but not virally suppressed'))]
g <- ggplot(subset(po)) + geom_bar(aes(fill=STAGE,y=M,x=FROM_BPLACE),stat='identity',position=position_dodge(width=0.9)) +
  geom_errorbar(aes(fill=STAGE,ymin=CL, ymax=CU,x=FROM_BPLACE),position=position_dodge(width=0.9), width=0.5, colour="black")	+
  scale_fill_npg('nrc') +
  #scale_fill_manual(name="Stage in cascade of\nlikely transmitter at\nputative infection time",values = c('Dutch-born'=pal[2],'Foreign-born'=pal[3])) +
  labs(x='Birthplace of recipient', y='Proportion of attributable\ninfections to place of birth', fill = 'Stage in cascade of likely transmitter at putative infection time') +
  theme_bw(base_size=28) +
  theme(legend.pos='bottom',
        #axis.text.y=element_blank(),
        strip.background=element_blank(),
        strip.text=element_blank(),
        axis.text.x = element_text(angle=45,vjust=0.7)) + #,
  guides(fill = guide_legend(nrow = 3,title.position='top'),by.col=T) +
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))
ggsave(file = paste0(outfile.base,'-flows_frombplace_stage.pdf'), g , w = 16, h = 12)
ggsave(file = paste0(outfile.base,'-flows_frombplace_stage.png'), g , w = 16, h = 12)


## plot ADJUSTED flows by birthplace and stage in cascade ----
cat(" \n --------------------------------  plot ADJUSTED flows by birthplace of case and stage of source -------------------------------- \n")

po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = 'tpair_prob' # should we be using tpair_prob now? unadjusted probs?
)
po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po <- data.table(po)
po[, PAIR_ID := as.integer(gsub(paste0('tpair_prob\\[([0-9]+)\\]'),'\\1',as.character(variable)))]
tmp <- subset(do, select = c('PAIR_ID','STAGE','TO_BPLACE','FROM_SEQUENCE_ID','TO_SEQUENCE_ID'))
po <- merge(po, tmp, by = 'PAIR_ID')
setnames(po,'value','rho')

# expected number of missing trsm pairs m_j(z) from region z to case j
# predict counts from negbin using r = number of obs individuals per recipient, p = sampling prob
#tmp2 <- subset(do,select=c('TO_SEQUENCE_ID','STAGE','TO_BPLACE'))
#tmp2 <- tmp2[, list(N=length(STAGE)),
#             by=c('TO_SEQUENCE_ID','STAGE','TO_BPLACE')]
tmp2 <- subset(do,select=c('TO_SEQUENCE_ID','FROM_SEQUENCE_ID','STAGE','TO_BPLACE','YEAR_OF_INF_EST'))
tmp2 <- tmp2[, list(N=length(FROM_SEQUENCE_ID)),
             by=c('TO_SEQUENCE_ID','STAGE','TO_BPLACE')]

#dmiss <- merge(subset(sp,select=c('LOC_BIRTH_POS','psi','psi_und')),tmp2,by.x='LOC_BIRTH_POS',by.y='TO_BPLACE')
dmiss <- cbind(tmp2,subset(spall,select=c('psi','psi_und')))
dmiss[STAGE=='undiag_first_year' | STAGE=='undiag_more_1_year', psi:= psi_und]

# calculate median w(z) and multiply by the number of missing sources in cat z for recipient j (per MC sample)
dw <- po[, list(rho_med_src = median(rho)), by = c('draw','STAGE')]
#dw <- merge(dmiss,dw,by.x='LOC_BIRTH_POS',by.y='STAGE',allow.cartesian = T) # adds rows per MC sample and per recipient missing obs
dw <- merge(dmiss,dw,by.x='STAGE',by.y='STAGE',allow.cartesian = T) # adds rows per MC sample and per recipient missing obs
# simulate missing sources per MC sample
#dw[, N_miss:= rnbinom(N,1,psi)]
dw <- dw[, list(STAGE=STAGE,
                TO_BPLACE=TO_BPLACE, psi=psi,
                TO_SEQUENCE_ID= TO_SEQUENCE_ID,
                N=N,
                draw=draw,
                rho_med_src=rho_med_src,
                N_miss= rnbinom(n=nrow(dw),size=N,prob=psi))]
#setnames(dw,'LOC_BIRTH_POS','STAGE')

po <- merge(po, dw, by=c('draw','TO_SEQUENCE_ID','STAGE'))

# calculate rho_kj = sum over observed rhos for incident case j
# sum over source cat z for unobserved sources for incident case j
tmp <- po[, list(rho_src = sum(rho),
                 mwz = sum(N_miss * rho_med_src)),by=c('draw','TO_SEQUENCE_ID')]
po <- merge(po,tmp, by=c('draw','TO_SEQUENCE_ID'))
po[rho_src + mwz==0, flag_zeroes:= 1]

# calculate new p_ij for obs pairs
# NB some values were NaN because rho, N_miss, rho_src and mwz are all 0. Set these to 0? or something small? 0.001?
dp_obs <- po[, list(p = rho/(rho_src + mwz)), by=c('draw','FROM_SEQUENCE_ID','TO_SEQUENCE_ID','STAGE')]
dp_obs[is.nan(p), p:=0] # <<< TODO: check this

# calculate new p_ij for missing pairs
dp_miss <- po[, list(p = rho_med_src / (rho_src + mwz),
                     N_miss = N_miss), by=c('draw','STAGE','TO_SEQUENCE_ID')]
dp_miss[is.nan(p), p:=0] # <<< TODO: check this
dp_miss[is.infinite(p), p:=0] # <<< some are inf because denominator is 0 but numerator is non-zero (bc N_miss=0)

# calculate p_j(x) = sum p_ij for obs pairs over obs sources i from bplace z + sum p_j for missing pairs over # missing sources from bplace z
tmp <- dp_obs[, list(p_ijz = sum(p)),by=c('draw','TO_SEQUENCE_ID','STAGE')]
tmp2 <- dp_miss[, list(p_jz = sum(p*N_miss)),by=c('draw','TO_SEQUENCE_ID','STAGE')]

tmp <- merge(tmp,tmp2,by=c('draw','STAGE','TO_SEQUENCE_ID'))
dp <- tmp[, list(pjx = p_ijz + p_jz), by=c('draw','STAGE','TO_SEQUENCE_ID')]

# adjust for sampling of incident cases
dp <- merge(dp,unique(subset(do,select=c('TO_SEQUENCE_ID','TO_BPLACE'))),by='TO_SEQUENCE_ID')
dp <- merge(dp, subset(sp,select=c('LOC_BIRTH_POS','psi')),by.x='TO_BPLACE',by.y='LOC_BIRTH_POS')
# calculate flows from source cat Z
po <- dp[, list(value = sum(pjx/psi)), by = c('draw','STAGE','TO_BPLACE')] # adjust recipients for sampling too
tmp <- po[, list(total = sum(value)), by = c('draw','TO_BPLACE')]
#tmp <- dp[, list(total = sum(pjx)), by = c('draw')]
po <- merge(po, tmp, by = c('draw','TO_BPLACE'))
po[, paf := value/total]
saveRDS(po,file=paste0(outfile.base,'-flows_fromstage_tobplace_adjusted_samplingbias_mcsamples','.RDS'))

po <- po[,
         list( q = quantile(paf, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('STAGE','TO_BPLACE')
]
po <- dcast.data.table(po, TO_BPLACE+STAGE~stat, value.var = 'q')
#po[, TO_BPLACE:= 'Overall']
saveRDS(po,file=paste0(outfile.base,'-flows_fromstage_tobplace_adjusted_samplingbias','.RDS'))


po[, STAGE:= factor(STAGE,
                          levels=c('undiag_first_year','undiag_more_1_year','diagnosed_unsuppressed'),
                          labels=c('Undiagnosed (1st year)',' Undiagnosed (>1 year)','Diagnosed but unsuppressed'))]
pal <- pal_npg("nrc")(4)[c(1,3,4)]

g1 <- ggplot(subset(po,TO_BPLACE=='Overall')) +
  geom_bar(aes(x=STAGE,y=M,fill=STAGE),stat='identity',position=position_dodge(width=0.9)) +
  geom_errorbar(aes(x=STAGE,ymin=CL, ymax=CU,fill=STAGE),position=position_dodge(width=0.9),width=0.5, colour="black")	+
  scale_fill_npg() +
  labs(x='Birthplace of\nlikely source',fill='Birthplace of likely\ntransmitter', y='Proportion of attributable infections\nto place of birth') +
  theme_bw(base_size=28) +
  theme(legend.pos='none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9)) #+ #,

ggsave(file = paste0(outfile.base,'-adjusted_flows_fromstage_toplace_imputemissingdata.pdf'),
       g1, w = 11, h = 8)
ggsave(file = paste0(outfile.base,'-adjusted_flows_fromstage_toplace_imputemissingdata.png'),
       g1, w = 11, h = 8)

po[, TO_BPLACE:= factor(TO_BPLACE,
                        levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                 'S. America &\n Caribbean','E. & C. Europe','MENA','Other'))]

g1 <- ggplot(subset(po)) +
  geom_bar(aes(x=TO_BPLACE,y=M,fill=STAGE),stat='identity',position=position_dodge(width=0.9)) +
  geom_errorbar(aes(x=TO_BPLACE,ymin=CL, ymax=CU,fill=STAGE),position=position_dodge(width=0.9),width=0.5, colour="black")	+
  scale_fill_npg() +
  labs(x='Birthplace of\nincident case',fill='Stage in cascade\nof likely transmitter', y='Proportion of attributable\ninfections to stage\nin cascade') +
  guides(fill = guide_legend(nrow = 3,title.position='top'),by.col=T) +
  theme_bw(base_size=18) +
  theme(legend.pos='right',
        #axis.title.x = element_blank(),
        axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9)) #+ #,

ggsave(file = paste0(outfile.base,'-adjusted_flows_fromstage_toplace_imputemissingdata.pdf'),
       g1, w = 9, h = 6)
ggsave(file = paste0(outfile.base,'-adjusted_flows_fromstage_toplace_imputemissingdata.png'),
       g1, w = 9, h = 6)

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
saveRDS(po,file=paste0(outfile.base,'-stratified_flows_fromstage_adjusted_samplingbias','.RDS'))


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

ggsave(file = paste0(outfile.base,'-stratified_flows_fromstage_adjusted_samplingbias.pdf'),
       g2, w = 9, h = 6)
ggsave(file = paste0(outfile.base,'-stratified_flows_fromstage_adjusted_samplingbias.png'),
       g2, w = 9, h = 6)



g <- ggarrange(g1 + rremove("xlab") + theme_bw()+ theme(axis.text.x = element_text(angle=45, vjust = 0.9,hjust=0.9),legend.pos='bottom'),
               g2+ theme(axis.text.x = element_text(angle=45, vjust = 0.9,hjust=0.9))+ rremove("xlab") + rremove("ylab"),
               ncol=2,widths=c(0.35,0.65),align='h',common.legend=T,legend = "bottom")
#g_bplace <- annotate_figure(g,bottom = text_grob("Birthplace of recipient",size=28))

ggsave(file = paste0(outfile.base,'-flows_fromstage_adjusted_samplingbias_panel.pdf'),
       g, w = 9, h = 6)
ggsave(file = paste0(outfile.base,'-flows_fromstage_adjusted_samplingbias_panel.png'),
       g, w = 9, h = 6)

