
## preamble ----
require(data.table)
require(lubridate)

if (1)
{
  args <- list(
    source_dir = '~/Documents/GitHub/transmission.flows.by.birthplace',
    indir = '~/Box\ Sync/Roadmap',
    analysis = 'analysis_220713',
    pairs.dir = '/Users/alexb/Documents/GitHub/source.attr.with.infection.time.fork/out_Amsterdam/agegps_updated_criteria_MSM-2010_2022',
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

source(file.path(args$source_dir, 'R', 'functions.R'))

cat(" \n --------------------------------  with arguments -------------------------------- \n")

indir.phsc <- file.path(args$indir, args$analysis, 'subgraphs')
infile.seq <-	file.path(args$indir, 'Data', 'data_220331/SHM_2201_ROADMAP_220331_tblLAB_seq.rda')
infile.bas <- file.path(args$indir, 'Data', 'data_220331','SHM_2201_ROADMAP_220331_tblBAS.csv')
infile.meta <- file.path(args$indir, args$analysis, 'misc', '220713_sequence_labels.rda')
infile.st <- file.path(args$indir, 'Data', 'data_220331/ROADMAP_220331_All_Taxa_withsubtype.rda')

## 2 columns: entire Amsterdam MSM cohort, sequenced MSM cohort, Amsterdam MSM with estimated infection date 2010-2021, 2010-2021 sequenced
## i.e.: first 2 cols relate to sources, second 2 relate to incident cases
# rows: Age, region of birth, year of diagnosis, HIV subtype (sequenced only), TSI?

## load sequence/subtype data ----
dsubgraphtaxa <- extract_subgraphs(indir.phsc)

load(infile.st)
st.a <- unique(subset(ds, select=c(FASTASampleCode,SUBTYPE)))
setnames(st.a,'FASTASampleCode','SEQ_LABEL')
st.a <- subset(st.a,grepl('Amst',SEQ_LABEL))

regex.tip.label <- '^([A-Za-z]+)_+([0-9]+)_([0-9-]+)_([A-Z0-9]+)'
st.a[, ID:= as.numeric(gsub(regex.tip.label,'\\2',SEQ_LABEL))]

load(infile.seq)
setkey(ds, PATIENT, SEQ_D)
first	<- ds[, list(SEQUENCE=SEQ_LABEL[1]), by='PATIENT']

st.a <- st.a[SEQ_LABEL %in% first$SEQUENCE]

## load meta data ----
dbas <- fread(infile.bas)
da <- subset(dbas,REG_EVR_AMSTERDAM==1 & MODE==1) # amsterdam MSM only
da[, SEQ:= PATIENT %in% dsubgraphtaxa$ID]
da <- merge(da,st.a,by.x='PATIENT',by.y='ID',all.x=T)
da[!is.na(SUBTYPE) & ! SUBTYPE %in% c('B','01_AE','02_AG','06_cpx','A1','C','D','G') , SUBTYPE:= 'Other']
da[, SUBTYPE:= factor(SUBTYPE,levels=c('B','01_AE','02_AG','A1','C','G','D','06_cpx','Other'),
                      labels=c('B','01AE','02AG','A1','C','G','D','06_cpx','Other'))]

## load infection dates ----
dinf <- data.table(read.csv(file.path('data_Ams',args$analysis,'Infection_date_est_rec.csv')))
setnames(dinf,c("id",'estsctodiagMedian','estsctodiagLL','estsctodiagUL'),c("PATIENT",'SER_TO_DIAG','SER_TO_DIAG_LL','SER_TO_DIAG_UL'))
dinf <- unique(dinf)
dinf[,DIAGNOSIS_DATE:= as.Date(dinf[,hiv_pos_d],format="%Y-%m-%d")]
dinf[,DIAGNOSIS_DATE_N:= hivc.db.Date2numeric(dinf[,DIAGNOSIS_DATE])]
dinf[,EST_INF_DATE:= DIAGNOSIS_DATE_N-SER_TO_DIAG]
dinf[,EST_INF_DATE:= format(date_decimal(EST_INF_DATE), "%Y-%m-%d")]
dinf[,YEAR_OF_INF_EST := year(EST_INF_DATE)]

# calculate median and 95\% CIs for Dutch/non-Dutch
dinf <- merge(dinf,subset(dbas,select=c('PATIENT','MODE','ORIGIN')),by='PATIENT',all.x=T)
dinf[, mwmb:= 'Foreign-born']
dinf[ORIGIN=='NL', mwmb:= 'Dutch-born']
med_tsi <- dinf[MODE==1 & YEAR_OF_INF_EST>=2010 & YEAR_OF_INF_EST<=2015, list( q = quantile(SER_TO_DIAG, probs = c(0.5, 0.025, 0.975),na.rm=T ),
                             stat = c('M','CL', 'CU')),
                by = c('mwmb')]
med_tsi <- dcast(med_tsi,mwmb~stat,value.var='q')
med_tsi[, L:= paste0(round(M,2),' [', round(CL,2),'-',round(CU,2),']')]
med_tsi[, L_months:= paste0(round(M*12,2),' [', round(CL*12,2),'-',round(CU*12,2),']')]
med_tsi

set(dinf,NULL,c('MODE','ORIGIN'),NULL)

## load bplace data ----
load(infile.seq)

load(infile.meta)

da <- merge(da,subset(dind,select=c('PATIENT','LOC_BIRTH')),by.x='PATIENT',by.y='PATIENT',all.x=T)

## merge data sets ----
da <- merge(da,dinf,by='PATIENT',all.x=T)
da[, AGE_AT_DIAG:= DIAGNOSIS_DATE_N - year(BIRTH_D)]
da[, AGE_AT_INF:= AGE_AT_DIAG - SER_TO_DIAG]
da[, AGE_GP_AT_DIAG:= cut(AGE_AT_DIAG, breaks=c(15,25,35,45,55,100),labels=c('15-24','25-34','35-44','45-54','55+'))]
da[, AGE_GP_AT_INF:= cut(AGE_AT_INF, breaks=c(15,25,35,45,55,100),labels=c('15-24','25-34','35-44','45-54','55+'))]
da[is.na(AGE_GP_AT_INF), AGE_GP_AT_INF:= 'Unknown']
da[, AGE_GP_AT_INF:= factor(AGE_GP_AT_INF,levels=c('15-24','25-34','35-44','45-54','55+','Unknown'))]

# categorise bplace
da[, BPLACE:="Other"]
da[LOC_BIRTH %in% c("WEurope","NorthAm","Oceania"), BPLACE:="W.Europe,\nN.America,Oceania"]
da[LOC_BIRTH %in% c("EEurope", "CEurope"), BPLACE:="E. & C. Europe"]
da[LOC_BIRTH %in% c("LaAmCar"), BPLACE:="S. America &\n Caribbean"]
da[LOC_BIRTH %in% c("DutchCarSuriname"), BPLACE:="Suriname &\nDutch Caribbean"]
da[LOC_BIRTH %in% c("MENA"), BPLACE:="MENA"]
da[ORIGIN=="NL", BPLACE:="Netherlands"]

da[, SEQ_LABEL:= NULL]
da <- subset(unique(da))


## make table ----
cat('--------- make table ----------')

tmp <- list()
tmp[['age']] <- da[, list(N_all=length(unique(na.omit(PATIENT))),
                          N_all_seq=length(unique(na.omit(PATIENT[!is.na(SUBTYPE)]))),
                          N_cases=length(unique(na.omit(PATIENT[YEAR_OF_INF_EST>=2010 & YEAR_OF_INF_EST<2022]))),
                          N_cases_seq=length(unique(na.omit(PATIENT[!is.na(SUBTYPE) & YEAR_OF_INF_EST>=2010 & YEAR_OF_INF_EST<2022])))),
                   by=c('AGE_GP_AT_INF')]

tmp[['age']] <- tmp[['age']][order(AGE_GP_AT_INF),]
setnames(tmp[['age']],'AGE_GP_AT_INF','GROUP')
tmp[['age']][, variable:= 'AGE_GP_AT_INF']

tmp[['bplace']] <- da[, list(N_all=length(na.omit(unique(PATIENT))),
                          N_all_seq=length(unique(na.omit(PATIENT[!is.na(SUBTYPE)]))),
                          N_cases=length(unique(na.omit(PATIENT[YEAR_OF_INF_EST>=2010 & YEAR_OF_INF_EST<2022]))),
                          N_cases_seq=length(unique(na.omit(PATIENT[!is.na(SUBTYPE) & YEAR_OF_INF_EST>=2010 & YEAR_OF_INF_EST<2022])))),
                   by=c('BPLACE')]

tmp[['bplace']] <- tmp[['bplace']][order(-N_all),]
setnames(tmp[['bplace']],'BPLACE','GROUP')
tmp[['bplace']][, variable:= 'BPLACE']


tmp[['subtype']] <- da[!is.na(SUBTYPE), list(N_all='-',
                             N_all_seq=length(unique(na.omit(PATIENT[!is.na(SUBTYPE)]))),
                             N_cases='-',
                             N_cases_seq=length(unique(na.omit(PATIENT[!is.na(SUBTYPE) & YEAR_OF_INF_EST>=2010 & YEAR_OF_INF_EST<2022])))),
                      by=c('SUBTYPE')]

tmp[['subtype']] <- tmp[['subtype']][order(SUBTYPE),]
setnames(tmp[['subtype']],'SUBTYPE','GROUP')
tmp[['subtype']][, variable:= 'SUBTYPE']

tmp[['total']] <- da[, list(variable='Total',
                            GROUP='-',
                            N_all=length(unique(na.omit(PATIENT))),
                                             N_all_seq=length(unique(na.omit(PATIENT[!is.na(SUBTYPE)]))),
                                             N_cases=length(unique(na.omit(PATIENT[YEAR_OF_INF_EST>=2010 & YEAR_OF_INF_EST<2022]))),
                                             N_cases_seq=length(unique(na.omit(PATIENT[!is.na(SUBTYPE) & YEAR_OF_INF_EST>=2010 & YEAR_OF_INF_EST<2022]))))]

tmp <- do.call(`rbind`,tmp)
tmp[variable!='SUBTYPE', N_all_seq_pct := paste0(N_all_seq, ' (', round(as.numeric(N_all_seq/as.integer(N_all))*100,0),'%)')]
tmp[variable!='SUBTYPE', N_cases_seq_pct := paste0(N_cases_seq, ' (', round(as.numeric(N_cases_seq/as.integer(N_cases))*100,0),'%)')]
tmp[is.na(N_all_seq_pct), N_all_seq_pct:= paste0(N_all_seq, ' (-)')]
tmp[is.na(N_cases_seq_pct), N_cases_seq_pct:= paste0(N_cases_seq, ' (-)')]
set(tmp,NULL,c('N_all_seq','N_cases_seq'),NULL)
tmp <- tmp[, c('variable','GROUP','N_all','N_all_seq_pct','N_cases','N_cases_seq_pct')]

saveRDS(tmp,file=paste0(outfile.base,'-table-characteristics.RDS'))

## table of birth place vs. subtype (2010-2021) ----

tmp <- da[, list(N_cases_seq=length(unique(na.omit(PATIENT[!is.na(SUBTYPE) & YEAR_OF_INF_EST>=2010 & YEAR_OF_INF_EST<2022])))),
                      by=c('BPLACE','SUBTYPE')]
tmp[, BPLACE:= factor(BPLACE,levels=c('Netherlands','W.Europe,\nN.America,Oceania',
                                      "S. America &\n Caribbean","Suriname &\nDutch Caribbean",
                                      "MENA" ,"E. & C. Europe" ,'Other'))]
tmp[, SUBTYPE:= factor(SUBTYPE,levels=c('B','01AE','02AG','A1','C','G','D','06_cpx','Other'),
                                               labels=c('B','01AE','02AG','A1','C','G','D','06cpx','Other'))]
tmp <- tmp[order(BPLACE,SUBTYPE),]

#tmp <- dcast(tmp,BPLACE~SUBTYPE,value.var='N_cases_seq')



g <- ggplot(subset(tmp,!is.na(SUBTYPE))) +
  geom_bar(aes(x=BPLACE,y=N_cases_seq,fill=SUBTYPE),stat="identity") +
  scale_fill_aaas() +
  labs(x='Birthplace of\nincident case',fill='Subtype', y='Number of sequenced\nincident cases (2010-2021)\namong Amsterdam MSM ') +
  theme_bw() +
  theme(legend.pos='right',
        axis.title.x = element_blank(),
        #axis.text.x = element_text(angle=40, vjust = 0.5)) + #,
        axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9)) #+ #,
  #coord_cartesian(ylim = c(0,1))
  #scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))
ggsave(file = paste0(outfile.base,'-bplace_subtype.pdf'), g, w = 6, h = 5)
ggsave(file = paste0(outfile.base,'-bplace_subtype.png'), g, w = 6, h = 5)


tmp[, ST_B:= 'Non-B']
tmp[SUBTYPE=='B', ST_B:= 'B']
tmp[, MG:= 'Foreign born']
tmp[BPLACE=='Netherlands', MG:= 'Dutch-born']
tmp <- tmp[, list(N=sum(N_cases_seq)),by=c('MG','ST_B')]
tmp <- dcast(tmp,MG~ST_B,value.var='N')
tmp[, B_pct:= B/(B+`Non-B`)*100]
tmp[, NonB_pct:= `Non-B`/(B+`Non-B`)*100]

saveRDS(tmp,file=paste0(outfile.base,'-bplace_subtype_cases.RDS'))

## plot sampling prob

spy <- readRDS(file=paste0(outfile.base,'-sampling_prob_byyear_cases','.RDS'))


## table of birth place vs. subtype (whole cohort) ----

tmp <- da[, list(N_cases_seq=length(unique(na.omit(PATIENT[!is.na(SUBTYPE)])))),
          by=c('BPLACE','SUBTYPE','YEAR_OF_INF_EST')]
tmp[, BPLACE:= factor(BPLACE,levels=c('Netherlands','W.Europe,\nN.America,Oceania',
                                      "S. America &\n Caribbean","Suriname &\nDutch Caribbean",
                                      "MENA" ,"E. & C. Europe" ,'Other'))]
tmp[, SUBTYPE:= factor(SUBTYPE,levels=c('B','01AE','02AG','A1','C','G','D','06_cpx','Other'),
                       labels=c('B','01AE','02AG','A1','C','G','D','06cpx','Other'))]
tmp <- tmp[order(BPLACE,YEAR_OF_INF_EST,SUBTYPE),]

#tmp <- dcast(tmp,BPLACE~SUBTYPE,value.var='N_cases_seq')


g <- ggplot(subset(tmp,!is.na(SUBTYPE))) +
  geom_bar(aes(x=YEAR_OF_INF_EST,y=N_cases_seq,fill=SUBTYPE),stat="identity",position="stack") +
  facet_wrap(BPLACE~.) +
  scale_fill_aaas() +
  labs(x='Birthplace of\nincident case',fill='Subtype', y='Number of sequenced\nincident cases (whole cohort)\namong Amsterdam MSM ') +
  theme_bw() +
  theme(legend.pos='right',
        axis.title.x = element_blank(),
        #axis.text.x = element_text(angle=40, vjust = 0.5)) + #,
        axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9)) #+ #,
#coord_cartesian(ylim = c(0,1))
#scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))
ggsave(file = paste0(outfile.base,'-bplace_subtype_year.pdf'), g, w = 6, h = 5)
ggsave(file = paste0(outfile.base,'-bplace_subtype_year.png'), g, w = 6, h = 5)

g <- ggplot(subset(tmp,!is.na(SUBTYPE))) +
  geom_bar(aes(x=YEAR_OF_INF_EST,y=N_cases_seq,fill=SUBTYPE),stat="identity",position="fill") +
  facet_wrap(BPLACE~.) +
  scale_fill_aaas() +
  labs(x='Birthplace of\nincident case',fill='Subtype', y='Number of sequenced\nincident cases (whole cohort)\namong Amsterdam MSM ') +
  theme_bw() +
  theme(legend.pos='right',
        axis.title.x = element_blank(),
        #axis.text.x = element_text(angle=40, vjust = 0.5)) + #,
        axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9)) #+ #,
#coord_cartesian(ylim = c(0,1))
#scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))
ggsave(file = paste0(outfile.base,'-bplace_subtype_year_prop.pdf'), g, w = 6, h = 5)
ggsave(file = paste0(outfile.base,'-bplace_subtype_year_prop.png'), g, w = 6, h = 5)


## sequenced (by subtype) and incident cases 1980-2021 ----

tmp <- da[, list(N_cases_seq=length(unique(na.omit(PATIENT[!is.na(SUBTYPE)])))),
          by=c('BPLACE','SUBTYPE','YEAR_OF_INF_EST')]
tmp[, BPLACE:= factor(BPLACE,levels=c('Netherlands','W.Europe,\nN.America,Oceania',
                                      "S. America &\n Caribbean","Suriname &\nDutch Caribbean",
                                      "MENA" ,"E. & C. Europe" ,'Other'))]
tmp[, SUBTYPE:= factor(SUBTYPE,levels=c('B','01AE','02AG','A1','C','G','D','06_cpx','Other'),
                       labels=c('B','01AE','02AG','A1','C','G','D','06cpx','Other'))]
tmp <- tmp[order(BPLACE,YEAR_OF_INF_EST,SUBTYPE),]
tmp <- subset(tmp,!is.na(SUBTYPE))

# total PLHIV
tmp2 <- da[!is.na(YEAR_OF_INF_EST), list(N_cases=length(unique(PATIENT))),
          by=c('BPLACE','YEAR_OF_INF_EST')]
tmp2[, BPLACE:= factor(BPLACE,levels=c('Netherlands','W.Europe,\nN.America,Oceania',
                                      "S. America &\n Caribbean","Suriname &\nDutch Caribbean",
                                      "MENA" ,"E. & C. Europe" ,'Other'))]
# correct for undiagnosed
ds <- readRDS(file=file.path(args$undiagnosed,paste0('p_undiagnosed_byyear_MC_samples_cohort_2010_2015_MSM.RDS')))
ds <- ds[,
         list( q = quantile(av_undiagnosed, probs = c(0.5, 0.025, 0.975) ),
               stat = c('M','CL', 'CU')
         ),
         by = c('year','mg')
]
ds <- dcast.data.table(ds, year+mg  ~stat, value.var = 'q')
setnames(ds,'M','av_undiagnosed')
ds <- subset(ds,select=c('year','mg','av_undiagnosed'))

dmap <- data.table(mwmb=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                          'S. America &\n Caribbean','E. & C. Europe','MENA','Other'),
                   mgid=c(1,2,3,4,5,6,7))
ds <- merge(ds,dmap,by.x='mg',by.y='mgid')
setnames(ds,'mwmb','migrant_group')
ds <- merge(tmp2,ds,by.x=c('BPLACE','YEAR_OF_INF_EST'),by.y=c('migrant_group','year'),all=T)
ds[is.na(av_undiagnosed), av_undiagnosed:=0] # set pre-2010 prob(undiagnosed) to 0
ds[is.na(N_cases), N_cases:=0] # set pre-2010 prob(undiagnosed) to 0
## calculate total infected
ds[, N_inf:= round(N_cases/(1-av_undiagnosed))]
ds[, cases:= 'Estimated incident\ncases']

g <- ggplot(subset(tmp)) +
  geom_bar(data=ds,aes(x=YEAR_OF_INF_EST,y=N_inf,col=cases),fill='white',size=0.3,stat="identity") +
  geom_bar(aes(x=YEAR_OF_INF_EST,y=N_cases_seq,fill=SUBTYPE),stat="identity",position="stack") +
  facet_wrap(BPLACE~.) +
  scale_fill_aaas() +
  scale_colour_manual(values='black') +
  labs(x='Birthplace of\nincident case',fill='Subtype', y='Number of sequenced\nincident cases (whole cohort)\namong Amsterdam MSM ',
       col='') +
  theme_bw() +
  theme(legend.pos='right',
        axis.title.x = element_blank(),
        #axis.text.x = element_text(angle=40, vjust = 0.5)) + #,
        axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9),
        strip.background=element_blank()) #+ #,
#coord_cartesian(ylim = c(0,1))
#scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))
ggsave(file = paste0(outfile.base,'-bplace_subtype_year_incident_cases.pdf'), g, w = 6.5, h = 5)
ggsave(file = paste0(outfile.base,'-bplace_subtype_year_incident_cases.png'), g, w = 6.5, h = 5)

#tmp <- merge(ds,subset(tmp,!is.na(SUBTYPE)),by=c('BPLACE','YEAR_OF_INF_EST'))

# sum the non-Bs
nonb <- tmp[, list(N_cases_nonb=sum(N_cases_seq[!is.na(SUBTYPE) & SUBTYPE!='B'])),
            by=c('YEAR_OF_INF_EST','BPLACE')]
nonb <- merge(ds,nonb,by=c('BPLACE','YEAR_OF_INF_EST'),all=T)
nonb <- subset(nonb,!is.na(YEAR_OF_INF_EST))
nonb[is.na(N_cases_nonb), N_cases_nonb:= 0]
nonb[, cases_b:= N_inf - N_cases_nonb]

tmp <- subset(tmp,!is.na(SUBTYPE))
tmp <- merge(tmp,nonb,by=c('BPLACE','YEAR_OF_INF_EST'),all=T)
tmp <- subset(tmp,!is.na(SUBTYPE) & !is.na(YEAR_OF_INF_EST))
tmp[SUBTYPE=='B', N_cases_seq:=cases_b]

g <- ggplot(subset(tmp,!is.na(SUBTYPE))) +
  geom_bar(aes(x=YEAR_OF_INF_EST,y=N_cases_seq,fill=SUBTYPE),stat="identity",position="fill") +
  facet_wrap(BPLACE~.) +
  scale_fill_aaas() +
  labs(x='Birthplace of\nincident case',fill='Subtype', y='Number of sequenced\nincident cases (whole cohort)\namong Amsterdam MSM ') +
  theme_bw() +
  theme(legend.pos='right',
        axis.title.x = element_blank(),
        #axis.text.x = element_text(angle=40, vjust = 0.5)) + #,
        axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9)) #+ #,
#coord_cartesian(ylim = c(0,1))
#scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))
ggsave(file = paste0(outfile.base,'-bplace_subtype_year_assume_missing_b_prop.pdf'), g, w = 6, h = 5)
ggsave(file = paste0(outfile.base,'-bplace_subtype_year_assume_missing_b_prop.png'), g, w = 6, h = 5)

## give numbers 2010-2021 by region ----

tmp <- da[, list(N_cases_seq=length(unique(na.omit(PATIENT[!is.na(SUBTYPE) & YEAR_OF_INF_EST>=2010 & YEAR_OF_INF_EST<2022])))),
          by=c('BPLACE','SUBTYPE')]
tmp[, BPLACE:= factor(BPLACE,levels=c('Netherlands','W.Europe,\nN.America,Oceania',
                                      "S. America &\n Caribbean","Suriname &\nDutch Caribbean",
                                      "MENA" ,"E. & C. Europe" ,'Other'))]
tmp[, SUBTYPE:= factor(SUBTYPE,levels=c('B','01AE','02AG','A1','C','G','D','06_cpx','Other'),
                       labels=c('B','01AE','02AG','A1','C','G','D','06cpx','Other'))]
tmp <- tmp[order(BPLACE,SUBTYPE),]

tmp <- subset(tmp,!is.na(SUBTYPE))
tmp[is.na(N_cases_seq), N_cases_seq:= 0]
tmp[, MG:= 'Foreign born']
tmp[BPLACE=='Netherlands', MG:= 'Netherlands']

tmp2 <- tmp[, list(N=sum(N_cases_seq)),by=c('MG','SUBTYPE')]
tmp3 <- tmp[, list(N_tot=sum(N_cases_seq)),by=c('MG')]
tmp2 <- merge(tmp2,tmp3,by='MG')
tmp2[, CL:= round(Hmisc::binconf(x=N,n=N_tot,return.df=T)[2]*N_tot,0)]
tmp2[, CU:= round(Hmisc::binconf(x=N,n=N_tot,return.df=T)[3]*N_tot,0)]
tmp2[, L:= paste0(N, ' [',CL,'-',CU,']')]
tmp2 <- dcast(tmp2,MG~SUBTYPE,value.var='L')

tmp3 <- tmp[, list(N_tot=sum(N_cases_seq)),by=c('BPLACE')]
tmp <- merge(tmp,tmp3,by='BPLACE')
tmp[, CL:= round(Hmisc::binconf(x=N_cases_seq,n=N_tot,return.df=T)[2]*N_tot,0)]
tmp[, CU:= round(Hmisc::binconf(x=N_cases_seq,n=N_tot,return.df=T)[3]*N_tot,0)]
tmp[, L:= paste0(N_cases_seq, ' [',CL,'-',CU,']')]
tmp <- dcast(tmp,BPLACE~SUBTYPE,value.var='L')

setnames(tmp2,'MG','BPLACE')
tmp <- rbind(tmp2,subset(tmp,BPLACE!='Netherlands'))
#tmp[is.na(tmp)] <- 0

#tmp[, total:= rowSums(tmp[,2:8])]


## proportions!

tmp <- da[, list(N_cases_seq=length(unique(na.omit(PATIENT[!is.na(SUBTYPE) & YEAR_OF_INF_EST>=2010 & YEAR_OF_INF_EST<2022])))),
          by=c('BPLACE','SUBTYPE')]
tmp[, BPLACE:= factor(BPLACE,levels=c('Netherlands','W.Europe,\nN.America,Oceania',
                                      "S. America &\n Caribbean","Suriname &\nDutch Caribbean",
                                      "MENA" ,"E. & C. Europe" ,'Other'))]
tmp[, SUBTYPE:= factor(SUBTYPE,levels=c('B','01AE','02AG','A1','C','G','D','06_cpx','Other'),
                       labels=c('B','01AE','02AG','A1','C','G','D','06cpx','Other'))]
tmp <- tmp[order(BPLACE,SUBTYPE),]

tmp <- subset(tmp,!is.na(SUBTYPE))
tmp[is.na(N_cases_seq), N_cases_seq:= 0]
tmp[, MG:= 'Foreign born']
tmp[BPLACE=='Netherlands', MG:= 'Netherlands']

tmp2 <- tmp[, list(N=sum(N_cases_seq)),by=c('MG','SUBTYPE')]
tmp3 <- tmp[, list(N_tot=sum(N_cases_seq)),by=c('MG')]
tmp2 <- merge(tmp2,tmp3,by='MG')
tmp2[, c('p','CL','CU'):= round(Hmisc::binconf(x=N,n=N_tot,return.df=T)*100,0)]
tmp2[, L:= paste0(p, '% [',CL,'-',CU,'%]')]
tmp2 <- dcast(tmp2,MG~SUBTYPE,value.var='L')

tmp3 <- tmp[, list(N_tot=sum(N_cases_seq)),by=c('BPLACE')]
tmp <- merge(tmp,tmp3,by='BPLACE')
tmp[, c('p','CL','CU'):= round(Hmisc::binconf(x=N_cases_seq,n=N_tot,return.df=T)*100,0)]
tmp[, L:= paste0(p, '% [',CL,'-',CU,'%]')]
tmp <- dcast(tmp,BPLACE~SUBTYPE,value.var='L')

setnames(tmp2,'MG','BPLACE')
tmp <- rbind(tmp2,subset(tmp,BPLACE!='Netherlands'))

## panel plot ----
### plot A: incident cases ----

tmp <- da[, list(N_cases_seq=length(unique(na.omit(PATIENT[!is.na(SUBTYPE)])))),
          by=c('BPLACE','SUBTYPE','YEAR_OF_INF_EST')]
tmp[, BPLACE:= factor(BPLACE,levels=c('Netherlands','W.Europe,\nN.America,Oceania',
                                      "S. America &\n Caribbean","Suriname &\nDutch Caribbean",
                                      "MENA" ,"E. & C. Europe" ,'Other'))]
tmp[, SUBTYPE:= factor(SUBTYPE,levels=c('B','01AE','02AG','A1','C','G','D','06_cpx','Other'),
                       labels=c('B','01AE','02AG','A1','C','G','D','06cpx','Other'))]
tmp <- tmp[order(BPLACE,YEAR_OF_INF_EST,SUBTYPE),]
tmp <- subset(tmp,!is.na(SUBTYPE))
# group migrants together
tmp[, MB:= 'Foreign-born']
tmp[BPLACE=='Netherlands', MB:= 'Dutch-born']
# group non-Bs together
tmp[, ST:= 'Non-B']
tmp[SUBTYPE=='B', ST:= 'B']

tmp <- tmp[, list(N_cases_seq=sum(N_cases_seq)),by=c('MB','ST','YEAR_OF_INF_EST')]


# total PLHIV
tmp2 <- da[!is.na(YEAR_OF_INF_EST), list(N_cases=length(unique(PATIENT))),
           by=c('BPLACE','YEAR_OF_INF_EST')]
tmp2[, BPLACE:= factor(BPLACE,levels=c('Netherlands','W.Europe,\nN.America,Oceania',
                                       "S. America &\n Caribbean","Suriname &\nDutch Caribbean",
                                       "MENA" ,"E. & C. Europe" ,'Other'))]
# correct for undiagnosed
ds <- readRDS(file=file.path(args$undiagnosed,paste0('p_undiagnosed_byyear_MC_samples_cohort_2010_2015_MSM.RDS')))
ds <- ds[,
         list( q = quantile(av_undiagnosed, probs = c(0.5, 0.025, 0.975) ),
               stat = c('M','CL', 'CU')
         ),
         by = c('year','mg')
]
ds <- dcast.data.table(ds, year+mg  ~stat, value.var = 'q')
setnames(ds,'M','av_undiagnosed')
ds <- subset(ds,select=c('year','mg','av_undiagnosed'))

dmap <- data.table(mwmb=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                          'S. America &\n Caribbean','E. & C. Europe','MENA','Other'),
                   mgid=c(1,2,3,4,5,6,7))
ds <- merge(ds,dmap,by.x='mg',by.y='mgid')
setnames(ds,'mwmb','migrant_group')
ds <- merge(tmp2,ds,by.x=c('BPLACE','YEAR_OF_INF_EST'),by.y=c('migrant_group','year'),all=T)
ds[is.na(av_undiagnosed), av_undiagnosed:=0] # set pre-2010 prob(undiagnosed) to 0
ds[is.na(N_cases), N_cases:=0] # set pre-2010 prob(undiagnosed) to 0
## calculate total infected
ds[, N_inf:= round(N_cases/(1-av_undiagnosed))]
ds[, cases:= 'Estimated incident\ncases']

# group migrants together
ds[, MB:= 'Foreign-born']
ds[BPLACE=='Netherlands', MB:= 'Dutch-born']

ds <- ds[, list(N_inf=sum(N_inf)),by=c('MB','YEAR_OF_INF_EST','cases')]

g_cases <- ggplot(subset(tmp)) +
  geom_bar(data=ds,aes(x=YEAR_OF_INF_EST,y=N_inf,col=cases),fill='white',size=0.3,stat="identity") +
  geom_bar(aes(x=YEAR_OF_INF_EST,y=N_cases_seq,fill=ST),stat="identity",position="stack") +
  facet_grid(MB~.) +
  scale_fill_aaas() +
  scale_colour_manual(values='black') +
  labs(x='Birthplace of\nincident case',fill='Subtype', y='Number of sequenced\nincident cases (whole cohort)\namong Amsterdam MSM ',
       col='') +
  theme_bw() +
  theme(legend.pos='right',
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9),
        strip.background=element_blank())

# colour by region of birth, shade subtype
# group the two datasets together
tmp2 <- dcast(tmp,MB+YEAR_OF_INF_EST~ST,value.var='N_cases_seq')
ds <- merge(ds,tmp2,by=c('MB','YEAR_OF_INF_EST'),all=T)
ds <- subset(ds,!is.na(YEAR_OF_INF_EST))
ds[is.na(B), B:= 0]
ds[is.na(`Non-B`), `Non-B`:= 0]
ds[, Unseq:= N_inf - B - `Non-B`]
set(ds,NULL,c('N_inf','cases'),NULL)
ds <- melt(ds,id.vars=c('MB','YEAR_OF_INF_EST'),variable.name='SUBTYPE')
ds[, SUBTYPE:= factor(SUBTYPE, levels=c('Unseq','B','Non-B'),labels=c('Unsequenced','B','Non-B'))]

#pal <- c('grey',pal_aaas('default')(2))
pal <- pal_npg('nrc')(4)[c(3,4)]

g_cases <- ggplot(subset(ds)) +
  geom_bar(aes(x=YEAR_OF_INF_EST,y=value,fill=MB,alpha=SUBTYPE),stat="identity",position="stack") +
  facet_grid(MB~.) +
  #scale_fill_aaas() +
  scale_fill_manual(values=pal) +
  scale_alpha_manual(values=c(0.4,0.6,1)) +
  scale_colour_manual(values='black') +
  labs(x='Birthplace of\nincident case',alpha='Subtype', fill='Birthplace',y='Number of incident cases\namong Amsterdam MSM ',
       col='') +
  guides(fill="none") +
  theme_bw() +
  theme(legend.pos='right',
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9),
        strip.background=element_blank())


### plot B: proportion of Non-Bs ----
# 1/ assuming all unsequenced were B
ds[, SUBTYPE2:=SUBTYPE]
ds[SUBTYPE=='Unsequenced', SUBTYPE2:='B']
ds_b <- ds[, list(allb=sum(value)),by=c('MB','YEAR_OF_INF_EST','SUBTYPE2')]

# 2/ assuming birthplace predicts B/Non-Bs (i.e. same proportion of B/Non-B)
ds[, YR_GP:= cut(YEAR_OF_INF_EST, breaks=c(1980,1985,1990,1995,2000,2005,2010,2015,
                                           2020,2025),
                                  labels=c('1980-1984','1985-1989','1990-1994','1995-1999','2000-2004',
                                           '2005-2009','2010-2014','2015-2019','2020-2024'),
                 include.lowest=T,right=F)]
tmp <- ds[SUBTYPE!='Unsequenced', list(p_nonB=sum(value[SUBTYPE=='Non-B']/sum(value))),
          by=c('MB' ,'YR_GP')]
ds <- merge(ds,tmp,by=c('MB' ,'YR_GP'))
ds2 <- dcast(ds,MB+YEAR_OF_INF_EST+p_nonB~SUBTYPE,value.var=c('value'))
ds2[, B:= round(B + (1-0.2841530)*Unsequenced,0)]
ds2[, `Non-B`:= round(`Non-B` + 0.2841530*Unsequenced,0)]
ds2 <- melt(ds2,id.vars=c('MB','YEAR_OF_INF_EST','p_nonB'),variable.name='SUBTYPE')
setnames(ds2,'value','prop_bplace')

setnames(ds_b,'SUBTYPE2','SUBTYPE')
ds <- merge(ds_b,subset(ds2,select=c('MB','YEAR_OF_INF_EST','SUBTYPE','prop_bplace')),by=c('MB','YEAR_OF_INF_EST','SUBTYPE'))
tmp <- melt(ds,id.vars=c('MB','YEAR_OF_INF_EST','SUBTYPE'),variable.name='assumption')
#tmp2 <- tmp[, list(SUBTYPE,
#                   p=value/sum(value)),by=c('MB','YEAR_OF_INF_EST','assumption')]
tmp2 <- tmp[, list(n_tot=sum(value)),by=c('MB','YEAR_OF_INF_EST','assumption')]
tmp <- merge(tmp,tmp2,by=c('MB','YEAR_OF_INF_EST','assumption'))
tmp[, c('p','CL','CU') := Hmisc::binconf(x=value,n=n_tot,return.df=T)]

tmp[, assumption:= factor(assumption,levels=c('prop_bplace','allb'),labels=c('Subtype predicted\nby birthplace','All unsequenced assumed\nB subtype'))]
#tmp2[is.nan(p), p]
# plot prop of non-Bs
pal <- pal_npg('nrc')(4)[c(3,4)]
g_props <- ggplot(subset(tmp,SUBTYPE=='Non-B')) + geom_point(aes(x=YEAR_OF_INF_EST,y=p,colour=MB),position=position_dodge(width=0.9)) +
  geom_errorbar(aes(x=YEAR_OF_INF_EST,ymin=CL, ymax=CU,fill=MB),position=position_dodge(width=0.9),width=0.5, colour="black")	+
  facet_grid(MB~assumption) +
  scale_colour_manual(values=pal) +
  labs(x='Year',fill='', y='Proportion of predicted\nnon-B cases',col='') +
  theme_bw() +
  theme(legend.pos='right',
        axis.title.x = element_blank(),
        #axis.text.x = element_text(angle=40, vjust = 0.5)) + #,
        axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9),
        strip.background=element_blank()) +#+ #,
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))


### plot C: time-shifting sources ----
po <- readRDS(file=paste0(outfile.base,'-rep_',replicate,'-adjusted_flows_mwmb_byyear_samplingofcases','.RDS'))

# for 3panel plot
g_srcs <- ggplot(subset(po,TO_BPLACE=='Overall')) +
  geom_bar(aes(x=YEAR_OF_INF_EST,y=M,fill=FROM_MIGRANT),position=position_dodge(width=0.9),stat='identity') +
  #geom_point(aes(x=YEAR_OF_INF_EST,y=M,colour=FROM_MIGRANT),position=position_dodge(width=0.9)) +
  geom_errorbar(aes(x=YEAR_OF_INF_EST,ymin=CL, ymax=CU,fill=FROM_MIGRANT),position=position_dodge(width=0.9),width=0.5, colour="black")	+
  scale_fill_manual(values=pal) +
  labs(x='',fill='Birthplace of likely\ntransmitter', y='Proportion of attributable infections\nto Dutch-born and foreign-born MSM',col='') +
  theme_bw() +
  theme(legend.pos='right',
        #axis.title.x = element_blank(),
        #axis.text.x = element_text(angle=40, vjust = 0.5)) + #,
        axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9)) +#+ #,
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))

# geom_point and facets
g_srcs <- ggplot(subset(po,TO_BPLACE=='Overall')) +
  #geom_bar(aes(x=YEAR_OF_INF_EST,y=M,fill=FROM_MIGRANT),position=position_dodge(width=0.9),stat='identity') +
  geom_point(aes(x=YEAR_OF_INF_EST,y=M,colour=FROM_MIGRANT),position=position_dodge(width=0.9)) +
  geom_errorbar(aes(x=YEAR_OF_INF_EST,ymin=CL, ymax=CU,fill=FROM_MIGRANT),position=position_dodge(width=0.9),width=0.5, colour="black")	+
  facet_grid(FROM_MIGRANT~.) +
  scale_fill_manual(values=pal) +
  scale_colour_manual(values=pal) +
  labs(x='',fill='Birthplace of likely\ntransmitter', y='Proportion of attributable infections\nto Dutch-born and foreign-born MSM',col='') +
  theme_bw() +
  theme(legend.pos='none',
        #axis.title.x = element_blank(),
        #axis.text.x = element_text(angle=40, vjust = 0.5)) + #,
        axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9),
        strip.background = element_blank()) +#+ #,
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))

# 2yr intervals
g_srcs <- ggplot(subset(po,TO_BPLACE=='Overall')) +
  #geom_bar(aes(x=YEAR_OF_INF_EST,y=M,fill=FROM_MIGRANT),position=position_dodge(width=0.9),stat='identity') +
  geom_point(aes(x=YEAR_GP,y=M,colour=FROM_MIGRANT),position=position_dodge(width=0.9)) +
  geom_errorbar(aes(x=YEAR_GP,ymin=CL, ymax=CU,fill=FROM_MIGRANT),position=position_dodge(width=0.9),width=0.5, colour="black")	+
  facet_grid(FROM_MIGRANT~.) +
  scale_fill_manual(values=pal) +
  scale_colour_manual(values=pal) +
  labs(x='',fill='Birthplace of likely\ntransmitter', y='Proportion of attributable infections\nto Dutch-born and foreign-born MSM',col='') +
  theme_bw() +
  theme(legend.pos='none',
        #axis.title.x = element_blank(),
        #axis.text.x = element_text(angle=40, vjust = 0.5)) + #,
        axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9),
        strip.background = element_blank()) +#+ #,
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))


g <- ggarrange(g_cases,g_props+theme(legend.position='none'),g_srcs,nrow=1,align='hv',
               labels='AUTO',font.label=list(size=14),widths=c(0.30,0.4,0.30))

ggsave(file=paste0(outfile.base,'-cases_subtypes_birthplace_sources_panel_v3.pdf'), g, w = 15, h = 6)
ggsave(file=paste0(outfile.base,'-cases_subtypes_birthplace_sources_panel_v3.png'), g, w = 15, h = 6)
