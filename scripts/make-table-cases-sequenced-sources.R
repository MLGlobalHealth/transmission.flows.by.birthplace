
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

tmp <- tmp[YEAR_OF_INF_EST %in% 2010:2021, list(N_cases_seq=sum(N_cases_seq)),by=c('BPLACE')]


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

tmp2 <- dcast(tmp,MB+YEAR_OF_INF_EST~ST,value.var='N_cases_seq')
ds <- merge(ds,tmp2,by=c('MB','YEAR_OF_INF_EST'),all=T)
ds <- subset(ds,!is.na(YEAR_OF_INF_EST))


# estimated number of incident cases per group
tab <- ds[YEAR_OF_INF_EST %in% 2010:2021, list(N_cases=sum(N_inf)),by='BPLACE']

# number phylogenetically observed incident cases in group
tmp

# number phylogenetically and epidemiologically possible transmission pairs identified by place of birth of incident case
do[,list(N_pairs=length(FROM_SEQUENCE_ID)),by='FROM_BPLACE']


# proportion of sources within each group
#-- from Dutch-born Amsterdam MSM
#-- etc

po <- readRDS(file=paste0(outfile.base,'-adjusted_flows_atob_samplingofcases_bplacecase_bplacesrc','.RDS'))

po <- dcast(po,FROM_BPLACE~TO_BPLACE,value.var='L')

