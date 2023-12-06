
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
require(lubridate)

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
  stopifnot(args_line[[7]]=='-analysis')
  stopifnot(args_line[[9]]=='-undiagnosed')
  stopifnot(args_line[[11]]=='-overwrite')

  args <- list()
  args[['source_dir']] <- args_line[[2]]
  args[['indir']] <- args_line[[4]]
  args[['outdir']] <- args_line[[6]]
  args[['analysis']] <- args_line[[8]]
  args[['undiagnosed']] <- args_line[[10]]
  args[['overwrite']] <- args_line[[12]]
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
## TODO FIX ARGS
args$analysis = args_dir$analysis
args$indir = args_dir$indir
#args$job_tag = args_dir$job_tag
args$undiagnosed = args_dir$undiagnosed
args$overwrite = args_dir$overwrite

tmp <- paste0(outfile.base,'-fitted_stan_model.rds')
cat("\n Read fitted dat ", tmp , "\n")
model_fit <- readRDS(file = tmp)

## load ethnicity data ----

load(infile.seq)

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

## calculate prevalence ----
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

### count the cases ----
# count the incident cases 2010-2021
# count the number who died each year (to remove in prevalence calculation)
tmp <- list()
tmp2 <- list()
bplaces <- unique(dinf$LOC_BIRTH_POS)
for(i in 1980:2021){
  tmp[[i]] <- list()
  tmp2[[i]] <- list()
  for(j in bplaces){
  tmp[[i]][[j]] <- dinf[,list(YEAR = i,
                         LOC_BIRTH_POS = j,
                         N_i=length(unique(na.omit(TO_SEQUENCE_ID[YEAR_OF_INF_EST== i & LOC_BIRTH_POS==j]))),
                         N_d=length(unique(na.omit(TO_SEQUENCE_ID[year(DEATH_D) == i & LOC_BIRTH_POS==j]))))]
  tmp2[[i]][[j]] <- dinf[,list(YEAR = i,
                          LOC_BIRTH_POS = j,
                          N_i_total=length(na.omit(TO_SEQUENCE_ID[YEAR_OF_INF_EST== i & LOC_BIRTH_POS==j])),
                          N_d_total=length(na.omit(TO_SEQUENCE_ID[year(DEATH_D) == i & LOC_BIRTH_POS==j])))]
  }
  tmp[[i]] <- do.call(`rbind`,tmp[[i]])
  tmp2[[i]] <- do.call(`rbind`,tmp2[[i]])
}
tmp <- do.call(`rbind`,tmp)
tmp2 <- do.call(`rbind`,tmp2)
tmp <- merge(tmp,tmp2,by=c('YEAR','LOC_BIRTH_POS'),all=T)

# plot new cases and deaths per year
ggplot(tmp) + geom_line(aes(x=YEAR,y=N_i_total),col='black') +
  geom_line(aes(x=YEAR,y=N_d_total),col='red') +
  labs(x='Year',y='Number of cases') +
  theme_bw()
if(args$overwrite){ggsave(file = paste0(outfile.base,'-incident_cases_deaths_byyear.png'), w = 5, h = 5)}

### adjust for undiagnosed ----
### load estimates of proportion undiagnosed

#ds <- readRDS(file=file.path(args$undiagnosed,paste0('p_undiagnosed_byyear_cohort_2010_2015_MSM.RDS')))
#ds <- dcast(ds,migrant_group+year~qlabel,value.var='p')
ds <- readRDS(file=file.path(args$undiagnosed,paste0('p_undiagnosed_byyear_MC_samples_cohort_2010_2015_MSM.RDS')))
dmap <- data.table(mwmb=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                          'S. America &\n Caribbean','E. & C. Europe','MENA','Other'),
                   mgid=c(1,2,3,4,5,6,7))
ds <- merge(ds,dmap,by.x='mg',by.y='mgid')
setnames(ds,'mwmb','migrant_group')
ds <- merge(tmp,ds,by.x=c('LOC_BIRTH_POS','YEAR'),by.y=c('migrant_group','year'),all=T)
ds[is.na(av_undiagnosed), av_undiagnosed:=0] # set pre-2010 prob(undiagnosed) to 0
## calculate total infected
ds[, N_inf:= round(N_i/(1-av_undiagnosed))]

## calculate prevalence 2010-2021
# sum the new cases acquired up to each year i and remove any cases which died up to year i
tmp <- list()
tmp2 <- list()
for(i in 1980:2021){
  tmp[[i]] <- ds[,list(YEAR = i,
                         N=sum(N_inf[YEAR<= i]) - sum(N_d[YEAR<=i])),
                   by=c('LOC_BIRTH_POS','iter')]
  tmp2[[i]] <- ds[,list(YEAR = i,
                          N_total=sum(N_inf[YEAR<= i]) - sum(N_d[YEAR<=i])),
  by='iter']
}
tmp <- do.call(`rbind`,tmp)
tmp2 <- do.call(`rbind`,tmp2)
tmp <- merge(tmp,tmp2,by=c('YEAR','iter'),all=T)
tmp[, prev:= N/N_total]
if(args$overwrite){saveRDS(tmp,file=paste0(outfile.base,'-PLHIV_byyear_MCsamples','.RDS'))}


### calculate weighted average ----
# generate weights of total PLHIV across the 11 years
wts <- tmp[LOC_BIRTH_POS=='Netherlands' & YEAR %in% 2010:2021, list(YEAR=YEAR, w=N_total/sum(N_total)),
           by='iter'] # only do for NL because N_total same for all birth regions
tmp <- merge(tmp,wts,by=c('YEAR','iter'))
dp <- tmp[, list(pct=sum(prev*w)),by=c('LOC_BIRTH_POS','iter')]
dp[, FROM_BPLACE:= factor(LOC_BIRTH_POS,
                          levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\n Caribbean','E. & C. Europe','MENA','Other'),
                          labels=c('Netherlands','W.Europe,\nN.America & Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\nNon-Dutch Caribbean','E. & C. Europe','MENA','Other'))]
if(args$overwrite){saveRDS(dp,file=paste0(outfile.base,'-prevalence_MC_samples','.RDS'))}

dp <- dp[,
         list( q = quantile(pct, probs = c(0.5, 0.025, 0.975) ),
               stat = c('M','CL', 'CU')
         ),
         by = c('FROM_BPLACE')
]
dp <- dcast.data.table(dp, FROM_BPLACE  ~stat, value.var = 'q')
if(args$overwrite){saveRDS(dp,file=paste0(outfile.base,'-prevalence_CIs','.RDS'))}

### plot prevalence ----
g_prev <- ggplot(subset(dp)) + geom_bar(aes(x=FROM_BPLACE,y=M,fill=FROM_BPLACE),stat='identity',position=position_dodge(width=0.9)) +
  #geom_errorbar(aes(x=FROM_BPLACE,ymin=CL, ymax=CU,fill=FROM_BPLACE),position=position_dodge(width=0.9), width=0.5, colour="black")	+
  scale_fill_npg() +
  labs(x='Birthplace of\nincident case',fill='Birthplace of likely\ntransmitter', y='Prevalence by place of birth \n(weighted average 2010-2021)') +
  theme_bw(base_size=28) +
  theme(legend.pos='none',
        axis.title.x = element_blank(),
        #axis.text.x = element_text(angle=40, vjust = 0.5)) + #,
        axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9)) +#+ #,
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))
if(args$overwrite){
  ggsave(file = paste0(outfile.base,'-rep_',replicate,'-contribution_to_prevalence.pdf'), g_prev, w = 23, h = 10)
  ggsave(file = paste0(outfile.base,'-rep_',replicate,'-contribution_to_prevalence.png'), g_prev, w = 16, h = 10)
}

### plot flows by birthplace ----
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
if(args$overwrite){saveRDS(po,file=paste0(outfile.base,'-flows_frombplace_adjusted_samplingbias_mcsamples','.RDS'))}

po <- po[,
         list( q = quantile(paf, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('FROM_BPLACE')
]
po <- dcast.data.table(po, FROM_BPLACE~stat, value.var = 'q')
setnames(po,'FROM_BPLACE','FROM_BPLACE')
po[, TO_BPLACE:= 'Overall']
if(args$overwrite){saveRDS(po,file=paste0(outfile.base,'-rep_',replicate,'-adjusted_flows_samplingofcases','.RDS'))}

po[, FROM_BPLACE:= factor(FROM_BPLACE,
                          levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\n Caribbean','E. & C. Europe','MENA','Other'),
                          labels=c('Netherlands','W.Europe,\nN.America & Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\nNon-Dutch Caribbean','E. & C. Europe','MENA','Other'))]

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
if(args$overwrite){
  ggsave(file = paste0(outfile.base,'-rep_',replicate,'-contribution_to_flows.pdf'), g1, w = 23, h = 10)
  ggsave(file = paste0(outfile.base,'-rep_',replicate,'-contribution_to_flows.png'), g1, w = 16, h = 10)
}

### plot flows over prevelance contributions ----
cat(" \n --------------------------------  plot flows over prevalence contributions -------------------------------- \n")

po <- readRDS(file=paste0(outfile.base,'-flows_frombplace_adjusted_samplingbias_mcsamples','.RDS'))

po[, FROM_BPLACE:= factor(FROM_BPLACE,
                          levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\n Caribbean','E. & C. Europe','MENA','Other'),
                          labels=c('Netherlands','W.Europe,\nN.America & Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\nNon-Dutch Caribbean','E. & C. Europe','MENA','Other'))]

po <- merge(po, dp, by='FROM_BPLACE')
#po[, contr_prev:= paf/pct]
po[, contr_prev:= paf/M]
po <- po[,
         list( q = quantile(contr_prev, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('FROM_BPLACE')
]
po <- dcast.data.table(po, FROM_BPLACE~stat, value.var = 'q')
po[, TO_BPLACE:= 'Overall']
if(args$overwrite){saveRDS(po,file=paste0(outfile.base,'-rep_',replicate,'-flows_frombplace_hivpos','.RDS'))}

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
if(args$overwrite){
  ggsave(file = paste0(outfile.base,'-rep_',replicate,'-contribution_to_flows_over_prev.pdf'), g2, w = 23, h = 10)
  ggsave(file = paste0(outfile.base,'-rep_',replicate,'-contribution_to_flows_over_prev.png'), g2, w = 16, h = 10)
}


## Make panel plot ----

leg <- get_legend(g2 + labs(fill='World region of birth') +   theme_bw(base_size=12) +
        theme(legend.position='bottom', plot.margin=unit(margin(-2,-2,-2,-2), "cm"),legend.margin=margin(c(0,0,0,0))) +
          guides(fill = guide_legend(nrow=3,title.position='top')))

g_l <- ggarrange(g1+theme_bw(base_size=11) + labs(y='\nContribution to transmissions') +
  theme(legend.pos='none',axis.title.x = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank()),
  labels=c('A'),font.label=list(size=14))
g_r <- ggarrange(ggarrange(
                  g_prev+theme_bw(base_size=11) + labs(y='\nContribution to prevalence\n2010-2021') + theme(legend.pos='none',
                                                                                                          axis.title.x = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank()),
                  g2+theme_bw(base_size=11) + labs(y='\nContribution to transmissions\n relative to contribution to prevalence\n2010-2021') + theme(legend.pos='none',
                                                                                                                                                  axis.title.x = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank()),
                ncol=2,nrow=1,align='v',labels=c('B','C'),font.label=list(size=14)),
                ggarrange(leg),ncol=1,nrow=2,align='hv',heights=c(0.5,0.5))
g <- ggarrange(g_l,g_r,ncol=2,align='hv',widths=c(0.35,0.65))

if(args$overwrite){
  ggsave(file = paste0(outfile.base,'-prev_flows_frombplace_contributions_panel_new_labs.pdf'),
       g, w = 10, h = 5)
  ggsave(file = paste0(outfile.base,'-prev_flows_frombplace_contributions_panel_new_labs.png'),
       g, w = 10, h = 5)
}
