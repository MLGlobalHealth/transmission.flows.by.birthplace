
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
require(lubridate)

if (1)
{
  args <- list(
    source_dir = '~/Documents/GitHub/transmission.flows.by.birthplace',
    #indir = '~/Box\ Sync/Roadmap/source_attribution',
    indir = '~/Box\ Sync/Roadmap',
    analysis = 'analysis_220713',
    pairs.dir = '/Users/alexb/Documents/GitHub/v/out_Amsterdam/agegps_updated_criteria_MSM-2010_2022',
    #outdir = '/Users/alexb/Documents/GitHub/source.attr.with.infection.time.fork/out_Amsterdam/mm_bgUnif_piGP_221027b-agegps_TE16_MSM-2010_2022-1665619',
    outdir = '/Users/alexb/Documents/GitHub/transmission.flows.by.birthplace/out_Amsterdam/mm_bgUnif_piGP_221027b-agegps_sensanalysis_210216_MSM-618873',
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
## TODO FIX ARGS
args$job_tag = 'agegps_sensanalysis_210216_MSM-2010_2022'
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
tmp <- subset(do, select = c('PAIR_ID','FROM_BPLACE','TO_BPLACE','FROM_MIGRANT','TO_MIGRANT','YEAR_OF_INF_EST'))
#setnames(tmp,'FROM_BPLACE','BPLACE')
po <- merge(po, tmp, by = 'PAIR_ID')
po <- merge(po, subset(spy,select=c('LOC_BIRTH_POS','YEAR_OF_INF_EST','psi')),
            by.x=c('TO_BPLACE','YEAR_OF_INF_EST'),by.y=c('LOC_BIRTH_POS','YEAR_OF_INF_EST'))
po <- po[, list(value = sum(value/psi)), by = c('draw','FROM_MIGRANT','YEAR_OF_INF_EST')]
tmp <- po[, list(total = sum(value)), by = c('draw','YEAR_OF_INF_EST')]
po <- merge(po, tmp, by = c('draw','YEAR_OF_INF_EST'))
po[, paf := value/total]
saveRDS(po,file=paste0(outfile.base,'-flows_frombplace_byyear_adjusted_samplingbias_mcsamples','.RDS'))

po <- po[,
         list( q = quantile(paf, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('FROM_MIGRANT','YEAR_OF_INF_EST')
]
po <- dcast.data.table(po, FROM_MIGRANT+YEAR_OF_INF_EST~stat, value.var = 'q')
po[, TO_BPLACE:= 'Overall']
saveRDS(po,file=paste0(outfile.base,'-rep_',replicate,'-adjusted_flows_mwmb_byyear_samplingofcases','.RDS'))

#po[, FROM_BPLACE:= factor(FROM_BPLACE,
#                          levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
#                                   'S. America &\n Caribbean','E. & C. Europe','MENA','Other'))]
pal <- pal_npg('nrc')(4)[c(3,4)]

# for 3panel plot
g1 <- ggplot(subset(po,TO_BPLACE=='Overall')) + geom_point(aes(x=YEAR_OF_INF_EST,y=M,colour=FROM_MIGRANT),position=position_dodge(width=0.9)) +
  geom_errorbar(aes(x=YEAR_OF_INF_EST,ymin=CL, ymax=CU,fill=FROM_MIGRANT),position=position_dodge(width=0.9),width=0.5, colour="black")	+
  #geom_smooth(aes(x=YEAR_OF_INF_EST,y=M,colour=FROM_MIGRANT),position=position_dodge(width=0.9)) +
  stat_smooth(aes(x=YEAR_OF_INF_EST,y=M,colour=FROM_MIGRANT),position=position_dodge(width=0.9),method = 'loess', span = 0.5,linewidth=0.7) +
  #scale_colour_npg() +
  scale_colour_manual(values=pal) +
  labs(x='Year',fill='Birthplace of likely\ntransmitter', y='Proportion of attributable infections\nto Dutch-born and foreign-born MSM',col='') +
  theme_bw() +
  theme(legend.pos='right',
        axis.title.x = element_blank(),
        #axis.text.x = element_text(angle=40, vjust = 0.5)) + #,
        axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9)) +#+ #,
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))
ggsave(file = paste0(outfile.base,'-contribution_to_flows_mwmb_byyear.pdf'), g1, w = 7, h = 4)
ggsave(file = paste0(outfile.base,'-contribution_to_flows_mwmb_byyear.png'), g1, w = 7, h = 4)


### recode into 3 year intervals? ----

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
tmp <- subset(do, select = c('PAIR_ID','FROM_BPLACE','TO_BPLACE','FROM_COUNTRY','TO_COUNTRY','FROM_SEQUENCE_ID','TO_SEQUENCE_ID','FROM_MIGRANT','TO_MIGRANT','YEAR_OF_INF_EST'))
#setnames(tmp,'FROM_BPLACE','BPLACE')
po <- merge(po, tmp, by = 'PAIR_ID')
po <- merge(po, subset(spy,select=c('LOC_BIRTH_POS','YEAR_OF_INF_EST','psi')),
            by.x=c('TO_BPLACE','YEAR_OF_INF_EST'),by.y=c('LOC_BIRTH_POS','YEAR_OF_INF_EST'))
po[, YEAR_GP:= cut(YEAR_OF_INF_EST, breaks=c(2010,2013,2016,2019,2021),labels=c('2010-2012','2013-2015','2016-2018','2019-2021'),include.lowest=T,right=F)]
po <- po[, list(value = sum(value/psi)), by = c('draw','FROM_MIGRANT','YEAR_GP')]
tmp <- po[, list(total = sum(value)), by = c('draw','YEAR_GP')]
po <- merge(po, tmp, by = c('draw','YEAR_GP'))
po[, paf := value/total]

po <- po[,
         list( q = quantile(paf, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('FROM_MIGRANT','YEAR_GP')
]
po <- dcast.data.table(po, FROM_MIGRANT+YEAR_GP~stat, value.var = 'q')
po[, TO_BPLACE:= 'Overall']


# for 3panel plot
g1 <- ggplot(subset(po,TO_BPLACE=='Overall')) + geom_point(aes(x=YEAR_GP,y=M,colour=FROM_MIGRANT),position=position_dodge(width=0.9)) +
  geom_errorbar(aes(x=YEAR_GP,ymin=CL, ymax=CU,fill=FROM_MIGRANT),position=position_dodge(width=0.9),width=0.5, colour="black")	+
  #geom_smooth(aes(x=YEAR_OF_INF_EST,y=M,colour=FROM_MIGRANT),position=position_dodge(width=0.9)) +
  stat_smooth(aes(x=YEAR_GP,y=M,colour=FROM_MIGRANT),position=position_dodge(width=0.9),method = 'loess', span = 0.5,linewidth=0.7) +
  scale_colour_npg() +
  labs(x='Year',fill='Birthplace of likely\ntransmitter', y='Proportion of attributable infections\nto Dutch-born and foreign-born MSM') +
  theme_bw() +
  theme(legend.pos='none',
        axis.title.x = element_blank(),
        #axis.text.x = element_text(angle=40, vjust = 0.5)) + #,
        axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9)) +#+ #,
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))
ggsave(file = paste0(outfile.base,'-contribution_to_flows_mwmb_by3years.pdf'), g1, w = 6, h = 4)
ggsave(file = paste0(outfile.base,'-contribution_to_flows_mwmb_by3years.png'), g1, w = 6, h = 4)

### recode into 2 year intervals? ----

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
tmp <- subset(do, select = c('PAIR_ID','FROM_BPLACE','TO_BPLACE','FROM_COUNTRY','TO_COUNTRY','FROM_SEQUENCE_ID','TO_SEQUENCE_ID','FROM_MIGRANT','TO_MIGRANT','YEAR_OF_INF_EST'))
#setnames(tmp,'FROM_BPLACE','BPLACE')
po <- merge(po, tmp, by = 'PAIR_ID')
po <- merge(po, subset(spy,select=c('LOC_BIRTH_POS','YEAR_OF_INF_EST','psi')),
            by.x=c('TO_BPLACE','YEAR_OF_INF_EST'),by.y=c('LOC_BIRTH_POS','YEAR_OF_INF_EST'))
po[, YEAR_GP:= cut(YEAR_OF_INF_EST, breaks=c(2010,2012,2014,2016,2018,2020,2022),
                   labels=c('2010-2011','2012-2013','2014-2015','2016-2017','2018-2019','2020-2021'),include.lowest=T,right=F)]
po <- po[, list(value = sum(value/psi)), by = c('draw','FROM_MIGRANT','YEAR_GP')]
tmp <- po[, list(total = sum(value)), by = c('draw','YEAR_GP')]
po <- merge(po, tmp, by = c('draw','YEAR_GP'))
po[, paf := value/total]

po <- po[,
         list( q = quantile(paf, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('FROM_MIGRANT','YEAR_GP')
]
po <- dcast.data.table(po, FROM_MIGRANT+YEAR_GP~stat, value.var = 'q')
po[, TO_BPLACE:= 'Overall']
saveRDS(po,file=paste0(outfile.base,'-adjusted_flows_mwmb_by2years_samplingofcases','.RDS'))


# for 3panel plot
g1 <- ggplot(subset(po,TO_BPLACE=='Overall')) + geom_point(aes(x=YEAR_GP,y=M,colour=FROM_MIGRANT),position=position_dodge(width=0.9)) +
  geom_errorbar(aes(x=YEAR_GP,ymin=CL, ymax=CU,fill=FROM_MIGRANT),position=position_dodge(width=0.9),width=0.5, colour="black")	+
  #geom_smooth(aes(x=YEAR_OF_INF_EST,y=M,colour=FROM_MIGRANT),position=position_dodge(width=0.9)) +
  stat_smooth(aes(x=YEAR_GP,y=M,colour=FROM_MIGRANT),position=position_dodge(width=0.9),method = 'loess', span = 0.5,linewidth=0.7) +
  scale_colour_npg() +
  labs(x='Year',fill='Birthplace of likely\ntransmitter', y='Proportion of attributable infections\nto Dutch-born and foreign-born MSM') +
  theme_bw() +
  theme(legend.pos='none',
        axis.title.x = element_blank(),
        #axis.text.x = element_text(angle=40, vjust = 0.5)) + #,
        axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9)) +#+ #,
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))
ggsave(file = paste0(outfile.base,'-contribution_to_flows_mwmb_by2years.pdf'), g1, w = 6, h = 4)
ggsave(file = paste0(outfile.base,'-contribution_to_flows_mwmb_by2years.png'), g1, w = 6, h = 4)


### recode into 5 year intervals? ----

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
tmp <- subset(do, select = c('PAIR_ID','FROM_BPLACE','TO_BPLACE','FROM_COUNTRY','TO_COUNTRY','FROM_SEQUENCE_ID','TO_SEQUENCE_ID','FROM_MIGRANT','TO_MIGRANT','YEAR_OF_INF_EST'))
#setnames(tmp,'FROM_BPLACE','BPLACE')
po <- merge(po, tmp, by = 'PAIR_ID')
po <- merge(po, subset(spy,select=c('LOC_BIRTH_POS','YEAR_OF_INF_EST','psi')),
            by.x=c('TO_BPLACE','YEAR_OF_INF_EST'),by.y=c('LOC_BIRTH_POS','YEAR_OF_INF_EST'))
#po[, YEAR_GP:= cut(YEAR_OF_INF_EST, breaks=c(2010,2016,2021),labels=c('2010-2015','2016-2021'),include.lowest=T,right=F)]
po[, YEAR_GP:= cut(YEAR_OF_INF_EST, breaks=c(2010,2014,2021),labels=c('2010-2013','2014-2021'),include.lowest=T,right=F)]
po <- po[, list(value = sum(value/psi)), by = c('draw','FROM_MIGRANT','YEAR_GP')]
tmp <- po[, list(total = sum(value)), by = c('draw','YEAR_GP')]
po <- merge(po, tmp, by = c('draw','YEAR_GP'))
po[, paf := value/total]

po <- po[,
         list( q = quantile(paf, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('FROM_MIGRANT','YEAR_GP')
]
po <- dcast.data.table(po, FROM_MIGRANT+YEAR_GP~stat, value.var = 'q')
po[, TO_MIGRANT:= 'All']

pal <- pal_npg("nrc")(4)[c(1,3,4)]

# for 3panel plot
g1 <- ggplot(subset(po,TO_MIGRANT=='All')) + geom_bar(aes(x=TO_MIGRANT,y=M,fill=FROM_MIGRANT),stat='identity',position=position_dodge(width=0.9)) +
  geom_errorbar(aes(x=TO_MIGRANT,ymin=CL, ymax=CU,fill=FROM_MIGRANT),position=position_dodge(width=0.9),width=0.5, colour="black")	+
  #geom_smooth(aes(x=YEAR_OF_INF_EST,y=M,colour=FROM_MIGRANT),position=position_dodge(width=0.9)) +
  #scale_colour_npg() +
  scale_fill_manual(name="Birthplace of\nlikely transmitter",values = c('Overall'='grey','Dutch-born'=pal[2],'Foreign-born'=pal[3])) +
  facet_wrap(YEAR_GP~.) +
  labs(x='Year',fill='Birthplace of likely\ntransmitter', y='Proportion of attributable infections\nto Dutch-born and foreign-born MSM') +
  theme_bw() +
  theme(legend.pos='none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9),
        #axis.text.x = element_blank(),
        strip.background=element_blank()) +#+ #,
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))
ggsave(file = paste0(outfile.base,'-contribution_to_flows_mwmb_by5years_2010_2014_2021.pdf'), g1, w = 6, h = 4)
ggsave(file = paste0(outfile.base,'-contribution_to_flows_mwmb_by5years_2010_2014_2021.png'), g1, w = 6, h = 4)


## stratified ----
po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = 'tpair_prob_w' #tpair_prob_w
)
po <- data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po <- data.table(po)
po[, PAIR_ID := as.integer(gsub(paste0('tpair_prob_w\\[([0-9]+)\\]'),'\\1',as.character(variable)))]
tmp <- subset(do, select = c('PAIR_ID','FROM_BPLACE','TO_BPLACE','FROM_COUNTRY','TO_COUNTRY','FROM_SEQUENCE_ID','TO_SEQUENCE_ID','FROM_MIGRANT','TO_MIGRANT','YEAR_OF_INF_EST'))
#setnames(tmp,'FROM_BPLACE','BPLACE')
po <- merge(po, tmp, by = 'PAIR_ID')
po <- merge(po, subset(spy,select=c('LOC_BIRTH_POS','YEAR_OF_INF_EST','psi')),
            by.x=c('TO_BPLACE','YEAR_OF_INF_EST'),by.y=c('LOC_BIRTH_POS','YEAR_OF_INF_EST'))
#po[, YEAR_GP:= cut(YEAR_OF_INF_EST, breaks=c(2010,2016,2021),labels=c('2010-2015','2016-2021'),include.lowest=T,right=F)]
po[, YEAR_GP:= cut(YEAR_OF_INF_EST, breaks=c(2010,2014,2021),labels=c('2010-2013','2014-2021'),include.lowest=T,right=F)]
po <- po[, list(value = sum(value/psi)), by = c('draw','FROM_MIGRANT','TO_MIGRANT','YEAR_GP')]
tmp <- po[, list(total = sum(value)), by = c('draw','YEAR_GP','TO_MIGRANT')]
po <- merge(po, tmp, by = c('draw','YEAR_GP','TO_MIGRANT'))
po[, paf := value/total]

po <- po[,
         list( q = quantile(paf, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('FROM_MIGRANT','TO_MIGRANT','YEAR_GP')
]
po <- dcast.data.table(po, FROM_MIGRANT+TO_MIGRANT+YEAR_GP~stat, value.var = 'q')
po[, TO_BPLACE:= 'Overall']

pal <- pal_npg("nrc")(4)[c(1,3,4)]

# for 3panel plot
g2 <- ggplot(subset(po,TO_BPLACE=='Overall')) + geom_bar(aes(x=TO_MIGRANT,y=M,fill=FROM_MIGRANT),stat='identity',position=position_dodge(width=0.9)) +
  geom_errorbar(aes(x=TO_MIGRANT,ymin=CL, ymax=CU,fill=FROM_MIGRANT),position=position_dodge(width=0.9),width=0.5, colour="black")	+
  #geom_smooth(aes(x=YEAR_OF_INF_EST,y=M,colour=FROM_MIGRANT),position=position_dodge(width=0.9)) +
  #scale_colour_npg() +
  scale_fill_manual(name="Birthplace of\nlikely transmitter",values = c('Overall'='grey','Dutch-born'=pal[2],'Foreign-born'=pal[3])) +
  facet_wrap(YEAR_GP~.) +
  labs(x='Birthplace of incident case',fill='Birthplace of likely\ntransmitter', y='Proportion of attributable infections\nto Dutch-born and foreign-born MSM') +
  theme_bw() +
  theme(legend.pos='bottom',
        #axis.title.x = element_blank(),
        #axis.text.x = element_text(angle=40, vjust = 0.5)) + #,
        axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9),
        strip.background=element_blank()) +#+ #,
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))
ggsave(file = paste0(outfile.base,'-stratified_contribution_to_flows_mwmb_by5years_2010_2014_2021.pdf'), g2, w = 6, h = 4)
ggsave(file = paste0(outfile.base,'-stratified_contribution_to_flows_mwmb_by5years_2010_2014_2021.png'), g2, w = 6, h = 4)



g <- ggarrange(g1,g2,ncol=1,align='hv',labels='AUTO',font.label=list(size=14),heights=c(0.41,0.59))
ggsave(file = paste0(outfile.base,'-contribution_to_flows_mwmb_by5years_2010_2014_2021_panel.pdf'), g, w = 6, h = 8)
ggsave(file = paste0(outfile.base,'-contribution_to_flows_mwmb_by5years_2010_2014_2021_panel.png'), g, w = 6, h = 8)

