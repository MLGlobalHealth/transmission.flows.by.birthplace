
## preamble ----
require(data.table)  # data mangling
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

## load bplace data ----
load(infile.seq)

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
ggsave(file = paste0(outfile.base,'-bplace_subtype_year.pdf'), g, w = 6, h = 5)
ggsave(file = paste0(outfile.base,'-bplace_subtype_year.png'), g, w = 6, h = 5)


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
