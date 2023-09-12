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

if (1)
{
  args <- list(
    source_dir = '~/Documents/GitHub/source.attr.with.infection.time.fork',
    #indir = '~/Box\ Sync/Roadmap/source_attribution',
    indir = '~/Box\ Sync/Roadmap',
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
args$analysis = 'analysis_220713'
args$indir = '~/Box\ Sync/Roadmap'

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

## plot case-adjusted flows by birthplace ----
cat(" \n --------------------------------  plot adjusted flows by birthplace -------------------------------- \n")
spy <- readRDS(file=paste0(outfile.base,'-sampling_prob_byyear_cases','.RDS'))


po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = 'tpair_prob_w'
)
po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po[, PAIR_ID := as.integer(gsub(paste0('tpair_prob_w\\[([0-9]+)\\]'),'\\1',as.character(variable)))]
tmp <- subset(do, select = c('PAIR_ID','FROM_MIGRANT','TO_BPLACE','TRANS_STAGE','YEAR_OF_INF_EST'))
po <- merge(po, tmp, by = 'PAIR_ID')
po[YEAR_OF_INF_EST<2016, PERIOD:= '2010-2015']
po[YEAR_OF_INF_EST>=2016, PERIOD:= '2016-2021']
po <- merge(po, subset(spy,select=c('LOC_BIRTH_POS','YEAR_OF_INF_EST','psi')),
            by.x=c('TO_BPLACE','YEAR_OF_INF_EST'),by.y=c('LOC_BIRTH_POS','YEAR_OF_INF_EST'))
po <- po[, list(value = sum(value/psi)), by = c('draw','FROM_MIGRANT','PERIOD')]
tmp <- po[, list(total = sum(value)), by = c('draw','PERIOD')]
po <- merge(po, tmp, by = c('draw','PERIOD'))
po[, paf := value/total]
po <- po[,
         list( q = quantile(paf, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('FROM_MIGRANT','PERIOD')
]
po <- dcast.data.table(po, PERIOD + FROM_MIGRANT~stat, value.var = 'q')
po[, TO_MIGRANT:= 'Overall']
saveRDS(po,file=paste0(outfile.base,'-rep_',replicate,'-PAF_all_pairs_period','.RDS'))

pal <- pal_npg("nrc")(4)[c(1,3,4)]

g1 <- ggplot(subset(po,TO_MIGRANT=='Overall')) + geom_bar(aes(fill=FROM_MIGRANT,y=M,x=TO_MIGRANT),stat='identity',position=position_dodge(width=0.9)) +
  geom_errorbar(aes(x=TO_MIGRANT,ymin=CL, ymax=CU,fill=FROM_MIGRANT),width=0.5, colour="black",position=position_dodge(width=0.9))	+
  scale_fill_manual(name="Birthplace of\nlikely transmitter",values = c('Overall'='grey','Dutch-born'=pal[2],'Foreign-born'=pal[3])) +
  facet_wrap(PERIOD~.) +
  labs(x='', y='Proportion of attributable\ninfections to place of birth') +
  theme_bw(base_size=28) +
  theme(legend.pos='none', strip.background=element_blank()) + #,
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))


# stratify by bplace of recipient
po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = 'tpair_prob_w'
)
po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po[, PAIR_ID := as.integer(gsub(paste0('tpair_prob_w\\[([0-9]+)\\]'),'\\1',as.character(variable)))]
tmp <- subset(do, select = c('PAIR_ID','FROM_MIGRANT','TO_MIGRANT','TO_BPLACE','TRANS_STAGE','YEAR_OF_INF_EST'))
po <- merge(po, tmp, by = 'PAIR_ID')
po[YEAR_OF_INF_EST<2016, PERIOD:= '2010-2015']
po[YEAR_OF_INF_EST>=2016, PERIOD:= '2016-2021']
po <- merge(po, subset(spy,select=c('LOC_BIRTH_POS','YEAR_OF_INF_EST','psi')),
            by.x=c('TO_BPLACE','YEAR_OF_INF_EST'),by.y=c('LOC_BIRTH_POS','YEAR_OF_INF_EST'))
po <- po[, list(value = sum(value/psi)), by = c('draw','FROM_MIGRANT','TO_MIGRANT','PERIOD')]
tmp <- po[, list(total = sum(value)), by = c('draw','TO_MIGRANT','PERIOD')]
po <- merge(po, tmp, by = c('draw','PERIOD','TO_MIGRANT'))
po[, paf := value/total]

po <- po[,
         list( q = quantile(paf, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('PERIOD','TO_MIGRANT','FROM_MIGRANT')
]
po <- dcast.data.table(po, PERIOD+TO_MIGRANT+FROM_MIGRANT~stat, value.var = 'q')
saveRDS(po,file=paste0(outfile.base,'-rep_',replicate,'-PAF_stratify_recipient_bplace_all_pairs_period','.RDS'))

g2 <- ggplot(subset(po,TO_MIGRANT!='Overall')) + geom_bar(aes(fill=FROM_MIGRANT,y=M,x=TO_MIGRANT),stat='identity',position=position_dodge(width=0.9)) +
  geom_errorbar(aes(fill=FROM_MIGRANT,ymin=CL, ymax=CU,x=TO_MIGRANT),position=position_dodge(width=0.9), width=0.5, colour="black")	+
  scale_fill_manual(name="Birthplace of\nlikely transmitter",values = c('Dutch-born'=pal[2],'Foreign-born'=pal[3])) +
  labs(x='Birthplace of recipient', y='Proportion of attributable\ninfections to place of birth') +
  facet_wrap(PERIOD~.) +
  theme_bw(base_size=28) +
  theme(legend.pos='bottom',
        #axis.text.y=element_blank(),
        strip.background=element_blank(),
        strip.text=element_blank()) + #,
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))

g <- ggarrange(g1 + rremove("ylab") , g2 + rremove("ylab"), ncol=1,heights=c(0.45,0.55),align='hv')
g_time <- annotate_figure(g, left = text_grob("Proportion of attributable infections\nto place of birth",size=28,rot = 90))

## combine plots ----

g_bplace <- g_bplace +  theme(plot.margin = margin(t=2,r=0,b=0,l=0, unit="cm"))
g_all <- ggarrange(g_bplace,g_time,ncol=1,heights=c(0.5,0.5),align='hv',
                   labels='AUTO',font.label=list(size=50),vjust=0.8)
ggsave(file = paste0(outfile.base,'-rep_',replicate,'-time_shifting_sources.pdf'), g_all, w = 16, h = 24)
ggsave(file = paste0(outfile.base,'-rep_',replicate,'-time_shifting_sources.png'), g_all, w = 16, h = 24)

