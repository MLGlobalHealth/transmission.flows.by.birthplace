require(data.table)
require(rstan)
require(bayesplot)
require(hexbin)
require(cmdstanr)
require(posterior)
require(truncdist)
require(scales)
require(ggsci)
require(tidyr)
require(reshape)
require(patchwork)
require(ggpubr)

## set up ----

args <- list(
  #source_dir= '/rds/general/user/ablenkin/home/git/bpm',
  #indir='/rds/general/project/ratmann_roadmap_data_analysis/live',
  #outdir= '/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model',
  #dir.ecdc='/rds/general/project/ratmann_roadmap_data_analysis/live/Data/Undiagnosed/ECDC_model',
  source_dir= '~/Documents/GitHub/transmission.flows.by.birthplace',
  indir='~/Box Sync/Roadmap',
  outdir= '~/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Roadmap/sources/ethnicity_analysis/undiagnosed',
  dir.ecdc='/Users/alexb/Box Sync/Roadmap/RQ1 Estimating introductions/undiagnosed/ECDC_model',
  stanModelFile= 'undiagnosed_211102',
  analysis= 'analysis_220713',
  #analysis= 'analysis_211101',
  #analysis= 'analysis_200917',
  hmc_stepsize= 0.02,
  hmc_num_samples= 15,
  hmc_num_warmup= 10,
  seed= 42,
  chain= 1,
  job_tag= 'dutch_v_migrant_2010_2015',
  sens=F,
  weights=F #'ECDC' unweighted for now since ECDC data only goes to 2019
)

## command line parsing if any
args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0)
{
  stopifnot(args_line[[1]]=='-source_dir')
  stopifnot(args_line[[3]]=='-stanModelFile')
  stopifnot(args_line[[5]]=='-analysis')
  stopifnot(args_line[[7]]=='-seed', !is.na(as.integer(args_line[[8]])))
  stopifnot(args_line[[9]]=='-indir')
  stopifnot(args_line[[11]]=='-outdir')
  stopifnot(args_line[[13]]=='-jobtag')
  stopifnot(args_line[[15]]=='-sens')
  stopifnot(args_line[[17]]=='-dir.ecdc')
  stopifnot(args_line[[19]]=='-weights')

  args <- list()
  args[['source_dir']] <- args_line[[2]]
  args[['stanModelFile']] <- args_line[[4]]
  args[['analysis']] <- args_line[[6]]
  args[['seed']] <- as.integer(args_line[[8]])
  args[['indir']] <- args_line[[10]]
  args[['outdir']] <- args_line[[12]]
  args[['job_tag']] <- args_line[[14]]
  args[['sens']] <- args_line[[16]]
  args[['dir.ecdc']] <- args_line[[18]]
  args[['weights']] <- args_line[[18]]
}

args$job_dir <- file.path(args$outdir,paste0(args$stanModelFile,'-',args$job_tag))
outpath <- file.path(args$outdir)
outdir <- file.path(args$outdir,paste0(args$stanModelFile,'-',args$job_tag))

dir.create( outdir )

args$trsm <- 'MSM'
job_tag <- args$job_tag

## make stan data ----
cat(" \n -------------------------------- \n Load data \n -------------------------------- \n")

# load data
file.seqlabels <- file.path(args$indir,args$analysis,'misc/220713_sequence_labels.rda')
infile.inftime <- file.path(args$indir,'Data/infection_time_estimates','220331_seq_update','roadmap_cd4_v3_est.csv')
if(args$sens!=F){
  infile.inftime <- file.path(args$indir,'Data/infection_time_estimates','220331_seq_update','roadmap_cd4_v3_est.csv')
}
geo.file <- file.path(args$indir,'misc/NEWGEO_220713.csv')

infile.ecdc.msm <- file.path(args$dir.ecdc,'AMS_MSM_LOCAL_MR_Result_main.csv')
infile.ecdc.nonmsm <- file.path(args$dir.ecdc,'AMS_NO_MSM_LOCAL_MR_Result_main.csv')

load(file.seqlabels)
dinf <- read.csv(infile.inftime,header=T)
geo <- data.table(read.csv(geo.file))

# load ECDC estimates
dw <- data.table(read.csv(infile.ecdc.msm))
dw[, trsm:='MSM']

cat(" \n -------------------------------- \n Generate weights \n -------------------------------- \n")

dw <- subset(dw, select=c(trsm,year,N_Inf_M), year %in% c(2014:2021))
dw <- dw[, list(year=year,w=N_Inf_M/sum(N_Inf_M)),by=c('trsm')]

cat(" \n -------------------------------- \n Define migrant groups \n -------------------------------- \n")

geo[geo$Alpha_2_code %in%c('AM','AZ','BY','GE','MD','RU','UA'),WRLD:='EEurope']
geo[geo$Alpha_2_code %in%c('AU','NZ'),WRLD:='Oceania']
setnames(geo,c('CNTRY','WRLD'),c('CNTRY_born','WRLD_born'))

dind <- data.table(dind)
dinf <- subset(dinf,select=c('id','estsctodiagMedian','hiv_pos_d'))
dinf <- unique(dinf)

dinf <- merge(dinf,subset(dind,select=c('PATIENT','TRANSM','BIRTH_CNTRY','HIV1_POS_D','INF_CNTRY','MIG_D')),by.x='id',by.y='PATIENT',all.x=T)
do <- data.table(dinf)
do[, time:=estsctodiagMedian]

do[, INF_D:=HIV1_POS_D - time]
do <- merge(do,subset(dind,select=c('PATIENT','ORIGIN')),by.x='id', by.y='PATIENT',all.x=T)
do <- merge(do,geo,by.x='ORIGIN',by.y='Alpha_2_code',all.x=T)

do[TRANSM=='MSM', mwmb:="Other"]
do[TRANSM=='MSM' & WRLD_born %in% c("WEurope","NorthAm","Oceania"), mwmb:="W.Europe,\nN.America,Oceania"]
do[TRANSM=='MSM' & WRLD_born %in% c("EEurope", "CEurope"), mwmb:="E. & C. Europe"]
do[TRANSM=='MSM' & WRLD_born %in% c("LaAmCar"), mwmb:="S. America &\n Caribbean"]
do[TRANSM=='MSM' & WRLD_born %in% c("DutchCarSuriname"), mwmb:="Suriname &\nDutch Caribbean"]
do[TRANSM=='MSM' & WRLD_born %in% c("MENA"), mwmb:="MENA"]
do[TRANSM=='MSM' & ORIGIN=="NL", mwmb:="Netherlands"]

# only keep MSM and those with a known birthplace
do <- subset(do, TRANSM=='MSM' & BIRTH_CNTRY!='Unknown')

# exclude individuals without an infection time estimate
do <- subset(do,!is.na(time))

# summarise number diagnosed to adjust for undiagnosed
do[, infdate:=INF_D]
do[is.na(INF_D) & HIV1_POS_D>=2014, infdate:=HIV1_POS_D] ## in model we use HIV pos date where there is no infection date - check this

do[, mwmb:= 'Migrant']
do[CNTRY_born=='Netherlands', mwmb:= 'Dutch']

n_diag <- do[, list(N_diag=length(unique(id[infdate>=2014]))),by=c('TRANSM','mwmb')]
n_diag <- n_diag[order(mwmb),]

# synthetic dataset includes an additional 3 years because we have diagnoses until end of 2021
da <- subset(do,do$INF_D>=2010 & do$INF_D<2016 & TRANSM %in% c('MSM'))

cat(" \n -------------------------------- \n Make summary table of characteristics \n -------------------------------- \n")

tab <- data.table(var1='TRANSM',var2='TOTAL',table(da$TRANSM))
tab[, V2:=V1]

tab <- rbind(tab,data.table(var1='TRANSM',var2='BIRTH_PLACE',table(da$mwmb,da$TRANSM)))
tab <- tab[, c('var1','var2','V1','V2','N')]

tab[var2=='TOTAL',V1:=var2]
tab <- subset(tab,!(N==0))
tab_all <- copy(tab)

# make synthetic subset of data
de <- subset(do,do$INF_D>=2010 & do$INF_D<2016 & TRANSM %in% c('HSX','MSM'))
dexcl <- da

tab <- data.table(var1='TRANSM',var2='TOTAL',table(dexcl$TRANSM))
tab[, V2:=V1]

tab <- rbind(tab,data.table(var1='TRANSM',var2='BIRTH_PLACE',table(dexcl$mwmb,dexcl$TRANSM)))
tab <- tab[, c('var1','var2','V1','V2','N')]
tab[var2=='TOTAL',V1:=var2]
tab <- subset(tab,!(N==0))
tab_excl <- copy(tab)

tab <- data.table(var1='TRANSM',var2='TOTAL',table(de$TRANSM))
tab[, V2:=V1]

tab <- rbind(tab,data.table(var1='TRANSM',var2='BIRTH_PLACE',table(de$mwmb,de$TRANSM)))
tab <- tab[, c('var1','var2','V1','V2','N')]

time <- de[, list(var='INFTIME',
                  V1='Estimated time to diagnosis (years)',
                  q=paste0(round(quantile(estsctodiagMedian,probs=0.5,na.rm=T),2),
                           " [",
                           round(quantile(estsctodiagMedian,probs=0.025,na.rm=T),2),
                           "-",
                           round(quantile(estsctodiagMedian,probs=0.975,na.rm=T),2),
                           "]")),
           by=c('TRANSM','mwmb')]
time_all <- de[, list(mwmb='TOTAL',
                      var='INFTIME',
                      V1='Estimated time to diagnosis (years)',
                      q=paste0(round(quantile(estsctodiagMedian,probs=0.5,na.rm=T),2),
                               " [",
                               round(quantile(estsctodiagMedian,probs=0.025,na.rm=T),2),
                               "-",
                               round(quantile(estsctodiagMedian,probs=0.975,na.rm=T),2),
                               "]")),
               by=c('TRANSM')]
time <- rbind(time,time_all)
time <- dcast(time,TRANSM+mwmb~V1,value.var='q')
tab[var2=='TOTAL',V1:=var2]
tab <- merge(tab,time,by.y=c('TRANSM','mwmb'),by.x=c('V2','V1'),all=T)
tab <- subset(tab,!(N==0 & is.na(`Estimated time to diagnosis (years)`)))
tab_incl <- copy(tab)

setnames(tab_all,'N','N_all')
setnames(tab_excl,'N','N_excluded')
setnames(tab_incl,'N','N_included')
tab <- merge(tab_all,tab_excl,by=c('var1','var2','V1','V2'),all=T)
tab <- merge(tab,tab_incl,by=c('var1','var2','V1','V2'),all=T)

tab[V1=='TOTAL', V1:='All']

tab$V2 <- factor(tab$V2,levels=c('MSM','HSX'),labels=c('Amsterdam MSM','Amsterdam heterosexual'))
tab <- tab[order(V2),]
#tab <- subset(tab,is.na(N_excluded)) # drop the pts infected outside NL
saveRDS(tab,file=file.path(outdir, paste0("characteristics_patients_undiagnosed.RDS")))
write.csv(tab,file=file.path(outdir, paste0("characteristics_patients_undiagnosed.csv")))

## do basic t test comparing TSIs ----

t.test(de$estsctodiagMedian[de$mwmb=='Dutch'],de$estsctodiagMedian[de$mwmb=='Migrant'],alternative = 'less',conf.level=0.95)
# 95%CI = -inf to -0.15 , p=0.001 ---> dutch-born have significantly shorter TSIs

## MSM model ----
cat(" \n -------------------------------- \n MSM model: Make stan data \n -------------------------------- \n")

args$trsm <- 'MSM'

dt <- subset(do,TRANSM==args$trsm & do$INF_D>=2010 & do$INF_D<2016)
N_diag <- subset(n_diag,TRANSM==args$trsm)

data_mg <- data.table(trsm='MSM',mgid=dt$mgid,migrant_group=dt$mwmb,time_to_diagnosis=dt$time)
saveRDS(data_mg,file=file.path(outdir, paste0('time_to_diagnosis_birthplace-',job_tag,"-",args$trsm,'.rds')))

q50 <- quantile(dt$time,probs=c(0.5))
q80 <- quantile(dt$time,probs=c(0.8))
stan_data <- list( n=nrow(dt), r=length(unique(dt$mwmb)), idx_to_obs_array=as.integer(factor(dt$mwmb)), y = dt$time,
                   log_wb_quantiles_priormean=c(log(q50),log(q80)-log(q50)),N_diag=N_diag$N_diag)
save(stan_data, file=file.path(outdir, paste0(job_tag,'-',args$trsm,'-stanin.RData')))

## init values
init_list <- list(
  list(wb_log_q50_overall=-1.2, wb_log_q80_q50_overall=0.7, wb_log_q50_sd=0.05,wb_log_q80_q50_sd=0.01),
  list(wb_log_q50_overall=-0.5, wb_log_q80_q50_overall=1.3, wb_log_q50_sd=0.2,wb_log_q80_q50_sd=0.1),
  list(wb_log_q50_overall=-0.1, wb_log_q80_q50_overall=1.6, wb_log_q50_sd=0.3,wb_log_q80_q50_sd=0.2)
)


cat(" \n -------------------------------- \n MSM model: run model \n -------------------------------- \n")

options(mc.cores=parallel::detectCores())
#options(mc.cores=1)
warmup <- 1000

# model using cmdstan
model = rstan::stan_model(file.path(args$source_dir,'stan_model_files','undiagnosed_211102.stan'))

fit = rstan::sampling(model,chains=3,data=stan_data,
                      warmup=500,iter=2000,
                      control=list(adapt_delta=.99),
                      init = init_list)
saveRDS(fit,file=file.path(outdir, paste0('stanfit_',job_tag,"_",args$trsm,'.rds')))

#	examine neff and rhat
fit.target.pars <- c('wb_log_q50_overall','wb_log_q50_grp','wb_log_q80_q50_overall','wb_log_q80_q50_grp','wb_log_q50_sd','wb_log_q80_q50_sd',
                     'wb_shape_grp[1]','wb_scale_grp[1]')
gqs.pars <- c('p_undiag_av','undiagnosed[1]')
summary <- rstan::monitor(rstan::extract(fit, pars=c(fit.target.pars,gqs.pars), permuted = FALSE, inc_warmup = TRUE))
print(summary,probs = c(0.025, 0.25, 0.5, 0.75, 0.975))

# model summary
mixture_su <- summary(fit)$summary
saveRDS(mixture_su,file=file.path(outdir, paste0('model_summary_fit_',job_tag,"_",args$trsm,'.rds')))

# get samples from the chains
samples <- rstan::extract(fit, inc_warmup = FALSE)
saveRDS(samples,file=file.path(outdir, paste0('samples_',job_tag,"_",args$trsm,'.rds')))

#	traces
color_scheme_set("mix-blue-red")
model.pars <- c('q50_overall','q50_r','q80_overall','q80_r','q50_sd','q80_sd')
p <- rstan::traceplot(fit, pars=c(fit.target.pars,'lp__'),inc_warmup=FALSE, ncol = 1)
pdf(file=file.path(outdir, paste0("undiagnosed_traces_sdprior_",job_tag,"_",args$trsm,".pdf")), w=10, h=20)
print(p)
dev.off()

# worst parameter
lp <- which(grepl("lp__",rownames(mixture_su)))
wrst <- rownames(mixture_su)[which.min(mixture_su[-lp,'n_eff'])]
small.neff <- mixture_su[wrst,]
write.csv(small.neff,file=paste0(outdir,'/undiag_untilmay2019',"_",args$trsm,'-pars-with-smallest-neff.csv'))
saveRDS(small.neff,file=paste0(outdir,'/undiag_untilmay2019',"_",args$trsm,'-pars-with-smallest-neff.rds'),version = 2)

p <- rstan::traceplot(fit, pars=wrst,inc_warmup=TRUE, ncol = 1)
pdf(file=file.path(outdir, paste0("undiagnosed_traces_sdprior-worst_par-",job_tag,"_",args$trsm,".pdf")), w=8, h=3)
print(p)
dev.off()

#	pair plots
model.pars <- c('q50_overall','q50_r[1]','q50_r[2]','q50_r[3]','q50_r[4]','q50_r[5]',
                'q80_overall','q80_r[1]','q80_r[2]','q80_r[3]','q80_r[4]','q80_r[5]',
                'q50_sd','q80_sd',"lp__")

p <- mcmc_pairs(rstan::extract(fit, pars=c(fit.target.pars,'lp__'), permuted=FALSE, inc_warmup=FALSE), diag_fun = "dens", off_diag_fun = "hex")
pdf(file=file.path(outdir, paste0("undiagnosed_mcmc_pairs_sdprior_",job_tag,"_",args$trsm,".pdf")), w=40, h=40)
print(p)
dev.off()

color_scheme_set("darkgray")
p <- mcmc_pairs(as.array(fit), np = nuts_params(fit), pars =c('wb_log_q50_overall','wb_log_q50_grp[1]','wb_log_q80_q50_overall',
                                                              'wb_log_q80_q50_grp[1]','wb_log_q50_sd','wb_log_q80_q50_sd',
                                                              'wb_shape_grp[1]','wb_scale_grp[1]',"lp__" ),
                diag_fun = "dens",off_diag_args = list(size = 0.75))
pdf(file=file.path(outdir, paste0("undiagnosed_mcmc_pairs_divergences_sdprior_",job_tag,"_",args$trsm,".pdf")), w=40, h=40)
print(p)
dev.off()


fit.target.pars <- c('wb_log_q50_overall','wb_log_q50_grp[1]','wb_log_q50_grp[2]','p_undiag_av[1]',
                     'p_undiag_av[2]','wb_log_q50_sd','wb_log_q80_q50_sd')
color_scheme_set("blue")
po <- rstan::extract(fit,inc_warmup=TRUE,permuted=FALSE)
p <- mcmc_intervals(po, pars = fit.target.pars)
p
ggsave(file=file.path(outdir,paste0('marginal_pairs_',job_tag,"_",args$trsm,'.pdf')), p, w=6, h=6)

## compare time since infection between groups ----

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

# calculate median
ds[, median:= scale*(log(2))^(1/shape)]
dm <- ds[, list(p=quantile(median,prob=c(0.025,0.5,0.975)),
                qlabel=c('p0.025','p0.5','p0.975')),by=c('trsm','mg')]
dm <- dcast(dm,trsm+mg~qlabel,value.var="p")
dm_all <- ds[, list(mg='All',
                    p=quantile(median,prob=c(0.025,0.5,0.975)),
                    qlabel=c('p0.025','p0.5','p0.975')),by=c('trsm')]
dm_all <- dcast(dm_all,trsm+mg~qlabel,value.var="p")
dm <- rbind(dm,dm_all)
saveRDS(dm,file=file.path(outdir,paste0('median_timetodiagnosis_',job_tag,"_",args$trsm,'.RDS')))

# test whether TSI in migrants > dutch-born
tmp <- dcast(ds,trsm+iter~mg,value.var='median')
tmp[, mig_longer_dutch:= ifelse(2>1,1,0)] # 100% of samples??



