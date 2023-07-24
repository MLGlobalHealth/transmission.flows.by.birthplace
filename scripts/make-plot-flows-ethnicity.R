
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

tmp <- list()
tmp2 <- list()
for(i in 2010:2021){
  tmp[[i]] <- dinf[,list(YEAR = i,
                         N=length(TO_SEQUENCE_ID[YEAR_OF_INF_EST<= i & (is.na(DEATH_D) | year(DEATH_D)>=YEAR_OF_INF_EST)])),
                   by=c('LOC_BIRTH_POS')]
  tmp2[[i]] <- dinf[,list(YEAR = i,
                          N_total=length(TO_SEQUENCE_ID[YEAR_OF_INF_EST<= i & (is.na(DEATH_D) | year(DEATH_D)>=YEAR_OF_INF_EST)])),
  ]
}
tmp <- do.call(`rbind`,tmp)
tmp2 <- do.call(`rbind`,tmp2)
tmp <- merge(tmp,tmp2,by=c('YEAR'),all=T)
tmp[, prev:= N/N_total]

## TODO! adjust for undiagnosed


# generate weights of total PLHIV across the 11 years
wts <- tmp[LOC_BIRTH_POS=='Netherlands', list(YEAR=YEAR, w=N_total/sum(N_total))] # only do for NL because N_total same for all birth regions
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

## plot flows by birthplace ----
cat(" \n --------------------------------  plot flows by birthplace -------------------------------- \n")

po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = 'tpair_prob_w'
)
po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po[, PAIR_ID := as.integer(gsub(paste0('tpair_prob_w\\[([0-9]+)\\]'),'\\1',as.character(variable)))]
tmp <- subset(do, select = c('PAIR_ID','FROM_BPLACE','TRANS_STAGE'))
setnames(tmp,'FROM_BPLACE','BPLACE')
po <- merge(po, tmp, by = 'PAIR_ID')
po <- po[, list(value = sum(value)), by = c('draw','BPLACE')]
tmp <- po[, list(total = sum(value)), by = c('draw')]
po <- merge(po, tmp, by = 'draw')
po[, paf := value/total]
po <- po[,
         list( q = quantile(paf, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('BPLACE')
]
po <- dcast.data.table(po, BPLACE~stat, value.var = 'q')
setnames(po,'BPLACE','FROM_BPLACE')
po[, TO_BPLACE:= 'Overall']
#saveRDS(po,file=paste0(outfile.base,'-rep_',replicate,'-PAF_bplace','.RDS'))

po <- readRDS(file=paste0(outfile.base,'-flows_frombplace_adjusted_samplingbias','.RDS'))

po[, FROM_BPLACE:= factor(FROM_BPLACE,
                          levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\n Caribbean','E. & C. Europe','MENA','Other'))]

pal <- pal_npg("nrc")(4)[c(1,3,4)]

g1 <- ggplot(subset(po,TO_BPLACE=='Overall')) + geom_bar(aes(x=TO_BPLACE,y=M,fill=FROM_BPLACE),stat='identity',position=position_dodge(width=0.9)) +
  geom_errorbar(aes(x=TO_BPLACE,ymin=CL, ymax=CU,fill=FROM_BPLACE),position=position_dodge(width=0.9),width=0.5, colour="black")	+
  scale_fill_npg() +
  labs(x='Birthplace of\nrecipient',fill='Birthplace of likely\ntransmitter', y='Proportion of attributable infections\nto place of birth') +
  theme_bw(base_size=28) +
  theme(legend.pos='none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=40, vjust = 0.5)) + #,
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))

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


#ggarrange(g1+theme(axis.text.x = element_text(angle=90, vjust = 0.5),
#                   plot.margin = unit(c(0,0,0,0), 'lines'),
#                   legend.pos='none') ,g2+ theme(axis.text.x = element_text(angle=90, vjust = 0.5),
#                                                                                 axis.text.y=element_blank(),
#                                                                                 axis.title = element_blank(),
#                                                                         plot.margin = unit(c(0,0,0,0), 'lines'),legend.pos='none'),
#          ncol=2,widths=c(0.35,0.65),align='hv')

#g_bplace <-ggarrange(g1+theme(axis.text.x = element_text(angle=90, vjust = 0.5),
#                   plot.margin = unit(c(0,1,0,0), 'lines')),
#          g2+ theme(axis.text.x = element_text(angle=90, vjust = 0.5),
#                    axis.text.y=element_blank(),
#                    axis.title.y = element_blank(),
#                    plot.margin = unit(c(0,0,1,0), 'lines')),
#          ncol=2,align='h',common.legend=T,legend = "bottom")
#g_bplace <- annotate_figure(g,bottom = text_grob("Birthplace of recipient",size=28))


## plot relative contribution ----

#infile.meta <- file.path(args$indir, analysis, 'misc', '220713_sequence_labels.rda')
#dind[, LOC_BIRTH_POS:="Other"]
#dind[LOC_BIRTH %in% c("WEurope","NorthAm","Oceania"), LOC_BIRTH_POS:="W.Europe,\nN.America,Oceania"]
#dind[LOC_BIRTH %in% c("EEurope", "CEurope"), LOC_BIRTH_POS:="E. & C. Europe"]
#dind[LOC_BIRTH %in% c("LaAmCar"), LOC_BIRTH_POS:="S. America &\n Caribbean"]
#dind[LOC_BIRTH %in% c("DutchCarSuriname"), LOC_BIRTH_POS:="Suriname &\nDutch Caribbean"]
#dind[LOC_BIRTH %in% c("MENA"), LOC_BIRTH_POS:="MENA"]
#dind[ORIGIN=="NL", LOC_BIRTH_POS:="Netherlands"]

#dp <- dind[CITY=='Amsterdam', list(pct=length(PATIENT)/nrow(dind[CITY=='Amsterdam'])),by='LOC_BIRTH_POS']

po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = 'tpair_prob_w'
)
po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po[, PAIR_ID := as.integer(gsub(paste0('tpair_prob_w\\[([0-9]+)\\]'),'\\1',as.character(variable)))]
tmp <- subset(do, select = c('PAIR_ID','FROM_BPLACE','TRANS_STAGE'))
setnames(tmp,'FROM_BPLACE','BPLACE')
po <- merge(po, tmp, by = 'PAIR_ID')
po <- po[, list(value = sum(value)), by = c('draw','BPLACE')]
tmp <- po[, list(total = sum(value)), by = c('draw')]
po <- merge(po, tmp, by = 'draw')
po[, paf := value/total]

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
  #coord_cartesian(ylim = c(0,1)) #+
  #scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))

## plot relative contribution among suppressed ----
# summarise birthplaces among all phylogenetically possible sources (i.e. all pairs)
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
new2 <- merge(new2,subset(dinf,select=c('TO_SEQUENCE_ID','LOC_BIRTH_POS')),by.x='PATIENT',by.y='TO_SEQUENCE_ID',all.x=T)
new2 <- subset(new2, !is.na(LOC_BIRTH_POS))

tmp <- new2[,list(N=length(PATIENT[pred_VL>200 & !is.na(pred_VL)])),
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

### plot unsuppressed ----
g_uns <- ggplot(subset(dp)) + geom_bar(aes(x=FROM_BPLACE,y=pct,fill=FROM_BPLACE),stat='identity',position=position_dodge(width=0.9)) +
  scale_fill_npg() +
  labs(x='Birthplace of\nincident case',fill='Birthplace of likely\ntransmitter', y='Proportion of unsuppressed\nby place of birth \n(weighted average 2010-2021)') +
  theme_bw(base_size=28) +
  theme(legend.pos='none',
        axis.title.x = element_blank(),
        #axis.text.x = element_text(angle=40, vjust = 0.5)) + #,
        axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9)) +#+ #,
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))

### plot ratio unsuppressed ----
po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = 'tpair_prob_w'
)
po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po[, PAIR_ID := as.integer(gsub(paste0('tpair_prob_w\\[([0-9]+)\\]'),'\\1',as.character(variable)))]
tmp <- subset(do, select = c('PAIR_ID','FROM_BPLACE','TRANS_STAGE'))
setnames(tmp,'FROM_BPLACE','BPLACE')
po <- merge(po, tmp, by = 'PAIR_ID')
po <- po[, list(value = sum(value)), by = c('draw','BPLACE')]
tmp <- po[, list(total = sum(value)), by = c('draw')]
po <- merge(po, tmp, by = 'draw')
po[, paf := value/total]

po <- readRDS(file=paste0(outfile.base,'-flows_frombplace_adjusted_samplingbias_mcsamples','.RDS'))

po <- merge(po, dp, by='FROM_BPLACE')
po[, contr_unsupp:= paf/pct]
po <- po[,
         list( q = quantile(contr_unsupp, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('FROM_BPLACE')
]
po <- dcast.data.table(po, FROM_BPLACE~stat, value.var = 'q')
po[, TO_BPLACE:= 'Overall']
po[, FROM_BPLACE:= factor(FROM_BPLACE,
                          levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\n Caribbean','E. & C. Europe','MENA','Other'))]
saveRDS(po,file=paste0(outfile.base,'-rep_',replicate,'-flows_frombplace_unsupp','.RDS'))

pal <- pal_npg("nrc")(4)[c(1,3,4)]

g3 <- ggplot(subset(po,TO_BPLACE=='Overall')) +
  geom_hline(yintercept=1,linetype=2) +
  geom_bar(aes(x=FROM_BPLACE,y=M,fill=FROM_BPLACE),stat='identity',position=position_dodge(width=0.9)) +
  geom_errorbar(aes(x=FROM_BPLACE,ymin=CL, ymax=CU,fill=FROM_BPLACE),position=position_dodge(width=0.9),width=0.5, colour="black")	+
  scale_fill_npg() +
  labs(x='Birthplace of\nincident case',fill='Birthplace of likely\ntransmitter', y='Contribution to transmission\nrelative to proportion unsuppressed\nphylogenetically possible sources') +
  theme_bw(base_size=28) +
  theme(legend.pos='none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9)) #+ #,

g <- ggarrange(g1+theme(axis.text.x = element_text(size=22,angle=60, vjust = 0.95,hjust = 0.9)),
               g2+theme(axis.text.x = element_text(size=22,angle=60, vjust = 0.95,hjust = 0.9)),
              g3+theme(axis.text.x = element_text(size=22,angle=60, vjust = 0.95,hjust = 0.9)),
              nrow=1, align='v',labels='AUTO',font.label=list(size=30),vjust=0.9)
g_bplace <- annotate_figure(g,bottom = text_grob("Birthplace of incident case",size=28))
ggsave(file = paste0(outfile.base,'-rep_',replicate,'-flows_frombplace_contributions.pdf'), g_bplace, w = 23, h = 10)
ggsave(file = paste0(outfile.base,'-rep_',replicate,'-flows_frombplace_contributions.png'), g_bplace, w = 16, h = 10)

g_top <- ggarrange(g_prev+theme_bw(base_size=11) + theme(legend.pos='none',
                                                     axis.title.x = element_blank(),axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9)),
                   g_uns+theme_bw(base_size=11) + theme(legend.pos='none',
                                                     axis.title.x = element_blank(),axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9)),
                   ncol=2,align='v',labels='AUTO',font.label=list(size=14),vjust=0.9)
g_bottom <- ggarrange(g1+theme_bw(base_size=11) + theme(legend.pos='none',
                                                 axis.title.x = element_blank(),axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9)),
               g2+theme_bw(base_size=11) + theme(legend.pos='none',
                                                 axis.title.x = element_blank(),axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9)),
               g3+theme_bw(base_size=11) + theme(legend.pos='none',
                                                 axis.title.x = element_blank(),axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9)),
               nrow=1, align='v',labels=c('C','D','E'),font.label=list(size=14),vjust=0.9)
g <- ggarrange(g_top,g_bottom,ncol=1,align='h')
g_bplace <- annotate_figure(g,bottom = text_grob("Birthplace of incident case",size=11))
ggsave(file = ('/Users/alexb/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Roadmap/sources/ethnicity_analysis/paper2_figures/prev_uns_flows_frombplace_contributions.pdf'),
       g_bplace, w = 11, h = 8)
ggsave(file = ('/Users/alexb/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Roadmap/sources/ethnicity_analysis/paper2_figures/prev_uns_flows_frombplace_contributions.png'),
       g_bplace, w = 11, h = 8)

## plot flows by birthplace and time period ----

po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = 'tpair_prob_w'
)
po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po[, PAIR_ID := as.integer(gsub(paste0('tpair_prob_w\\[([0-9]+)\\]'),'\\1',as.character(variable)))]
tmp <- subset(do, select = c('PAIR_ID','FROM_MIGRANT','TRANS_STAGE','YEAR_OF_INF_EST'))
po <- merge(po, tmp, by = 'PAIR_ID')
po[YEAR_OF_INF_EST<2016, PERIOD:= '2010-2015']
po[YEAR_OF_INF_EST>=2016, PERIOD:= '2016-2021']
po <- po[, list(value = sum(value)), by = c('draw','FROM_MIGRANT','PERIOD')]
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
tmp <- subset(do, select = c('PAIR_ID','FROM_MIGRANT','TO_MIGRANT','TRANS_STAGE','YEAR_OF_INF_EST'))
po <- merge(po, tmp, by = 'PAIR_ID')
po[YEAR_OF_INF_EST<2016, PERIOD:= '2010-2015']
po[YEAR_OF_INF_EST>=2016, PERIOD:= '2016-2021']
po <- po[, list(value = sum(value)), by = c('draw','FROM_MIGRANT','TO_MIGRANT','PERIOD')]
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
ggsave(file = paste0(outfile.base,'-rep_',replicate,'-PAF_panel_AB_lab.pdf'), g_all, w = 16, h = 24)
ggsave(file = paste0(outfile.base,'-rep_',replicate,'-PAF_panel_AB_lab.png'), g_all, w = 16, h = 24)

