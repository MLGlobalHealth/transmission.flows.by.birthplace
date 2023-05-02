
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
sp[, psi:= N_seq/N]

# calculate using birthplaces among incident cases in formulated pairs
#tmp2 <- unique(subset(do,select=c('TO_SEQUENCE_ID','TO_BPLACE')))
#tmp2 <- tmp2[, list(N_samp=length(unique(TO_SEQUENCE_ID)),
#                  pct_samp=length(unique(TO_SEQUENCE_ID))/nrow(tmp2)),
#           by='TO_BPLACE']

#sp <- merge(tmp,tmp2,by.x='LOC_BIRTH_POS',by.y='TO_BPLACE')
#sp[, psi:= N_samp/N]


## plot adjusted flows by birthplace ----
cat(" \n --------------------------------  plot adjusted flows by birthplace -------------------------------- \n")

po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = 'tpair_prob_w'
)
po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po[, PAIR_ID := as.integer(gsub(paste0('tpair_prob_w\\[([0-9]+)\\]'),'\\1',as.character(variable)))]
tmp <- subset(do, select = c('PAIR_ID','FROM_BPLACE','TO_BPLACE','TRANS_STAGE'))
#setnames(tmp,'FROM_BPLACE','BPLACE')
po <- merge(po, tmp, by = 'PAIR_ID')
po <- merge(po, subset(sp,select=c('LOC_BIRTH_POS','psi')),by.x='TO_BPLACE',by.y='LOC_BIRTH_POS')
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
saveRDS(po,file=paste0(outfile.base,'-rep_',replicate,'-adjusted_flows_tobplace','.RDS'))

po[, FROM_BPLACE:= factor(FROM_BPLACE,
                          levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\n Caribbean','E. & C. Europe','MENA','Other'))]

g <- ggplot(subset(po,TO_BPLACE=='Overall')) + geom_bar(aes(x=FROM_BPLACE,y=M,fill=FROM_BPLACE),stat='identity',position=position_dodge(width=0.9)) +
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

ggsave(file = paste0(outfile.base,'-adjusted_flows_tobplace_contributions.pdf'),
       g, w = 11, h = 8)
ggsave(file = paste0(outfile.base,'-adjusted_flows_tobplace_contributions.png'),
       g, w = 11, h = 8)


## adjust for sampling bias of sources ----

po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = 'tpair_prob_w' # should we be using tpair_prob now? unadjusted probs?
)
po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po[, PAIR_ID := as.integer(gsub(paste0('tpair_prob_w\\[([0-9]+)\\]'),'\\1',as.character(variable)))]
tmp <- subset(do, select = c('PAIR_ID','FROM_BPLACE','TO_BPLACE','FROM_SEQUENCE_ID','TO_SEQUENCE_ID'))
po <- merge(po, tmp, by = 'PAIR_ID')
setnames(po,'value','rho')

# expected number of missing trsm pairs m_j(z) from region z to case j
# predict counts from negbin using r = number of obs individuals per recipient, p = sampling prob
tmp2 <- subset(do,select=c('TO_SEQUENCE_ID','FROM_BPLACE'))
tmp2 <- tmp2[, list(N=length(FROM_BPLACE)),
             by=c('TO_SEQUENCE_ID','FROM_BPLACE')]

dmiss <- merge(subset(sp,select=c('LOC_BIRTH_POS','psi')),tmp2,by.x='LOC_BIRTH_POS',by.y='FROM_BPLACE')

# calculate median w(z) and multiply by the number of missing sources in cat z for recipient j (per MC sample)
dw <- po[, list(rho_med_src = median(rho)), by = c('draw','FROM_BPLACE')]
dw <- merge(dmiss,dw,by.x='LOC_BIRTH_POS',by.y='FROM_BPLACE',allow.cartesian = T) # adds rows per MC sample and per recipient missing obs
# simulate missing sources per MC sample
dw[, N_miss:= rnbinom(N,1,psi)]
setnames(dw,'LOC_BIRTH_POS','FROM_BPLACE')

po <- merge(po, dw, by=c('draw','TO_SEQUENCE_ID','FROM_BPLACE'))

# calculate rho_kj = sum over observed rhos for incident case j
# sum over source cat z for unobserved sources for incident case j
tmp <- po[, list(rho_src = sum(rho),
                 mwz = sum(N_miss * rho_med_src)),by=c('draw','TO_SEQUENCE_ID')]
po <- merge(po,tmp, by=c('draw','TO_SEQUENCE_ID'))
po[rho_src + mwz==0, flag_zeroes:= 1]

# calculate new p_ij for obs pairs
# NB some values were NaN because rho, N_miss, rho_src and mwz are all 0. Set these to 0? or something small? 0.001?
dp_obs <- po[, list(p = rho/(rho_src + mwz)), by=c('draw','FROM_SEQUENCE_ID','TO_SEQUENCE_ID','FROM_BPLACE')]
dp_obs[is.nan(p), p:=0] # <<< TODO: check this

# calculate new p_ij for missing pairs
dp_miss <- po[, list(p = rho_med_src / (rho_src + mwz),
                     N_miss = N_miss), by=c('draw','FROM_BPLACE','TO_SEQUENCE_ID')]
dp_miss[is.nan(p), p:=0] # <<< TODO: check this
dp_miss[is.infinite(p), p:=0] # <<< some are inf because denominator is 0 but numerator is non-zero (bc N_miss=0)

# calculate p_j(x) = sum p_ij for obs pairs over obs sources i from bplace z + sum p_j for missing pairs over # missing sources from bplace z
tmp <- dp_obs[, list(p_ijz = sum(p)),by=c('draw','TO_SEQUENCE_ID','FROM_BPLACE')]
tmp2 <- dp_miss[, list(p_jz = sum(p*N_miss)),by=c('draw','TO_SEQUENCE_ID','FROM_BPLACE')]

tmp <- merge(tmp,tmp2,by=c('draw','FROM_BPLACE','TO_SEQUENCE_ID'))
dp <- tmp[, list(pjx = p_ijz + p_jz), by=c('draw','FROM_BPLACE','TO_SEQUENCE_ID')]

# adjust for sampling of incident cases
dp <- merge(dp,unique(subset(do,select=c('TO_SEQUENCE_ID','TO_BPLACE'))),by='TO_SEQUENCE_ID')
dp <- merge(dp, subset(sp,select=c('LOC_BIRTH_POS','psi')),by.x='TO_BPLACE',by.y='LOC_BIRTH_POS')
# calculate flows from source cat Z
po <- dp[, list(value = sum(pjx/psi)), by = c('draw','FROM_BPLACE')]
tmp <- po[, list(total = sum(value)), by = c('draw')]
#tmp <- dp[, list(total = sum(pjx)), by = c('draw')]
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

