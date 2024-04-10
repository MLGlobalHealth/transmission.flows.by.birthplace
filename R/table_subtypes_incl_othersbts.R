
require(plyr)
require(dplyr)
require(tidyr)
require(ape)
require(data.table)

dir.name		<- '~/Box Sync/Roadmap/Data/data_220331'
outfile.matches      <- file.path(dir.name,'ROADMAP_220331_All_Taxa_withsubtype.rda')
outfile.trimseqs.pruned <- file.path(dir.name,'ROADMAP_220331_All_Sequences_LANL_aligned_curated.fasta')
infile.bas <- file.path('/Users/alexb/Box Sync/Roadmap', 'Data', 'data_220331','SHM_2201_ROADMAP_220331_tblBAS.csv')
infile.georeg <- file.path('~/Box Sync/Roadmap','misc','NEWGEO_220713.csv')
dgeo <- read.csv(infile.georeg,header=T)
dgeo <- dgeo %>% select(c('Alpha_2_code', 'Alpha_3_code','CNTRY','WRLD')) %>%
  dplyr::rename('Alpha.2.code'='Alpha_2_code') %>% dplyr::rename('Alpha.3.code'='Alpha_3_code') %>%
  mutate(Alpha.2.code:=as.character(Alpha.2.code),
         WRLD:= as.character(WRLD),
         CNTRY:= as.character(CNTRY))

dbas <- data.table(read.csv(infile.bas))

seqs.all <- read.dna(outfile.trimseqs.pruned,format='fa')
load(outfile.matches)

# remove some taxa which had very long branches in trees
#if(0){ # skip this! want to summarise all LANL selected, ignoring poor alignments
# only remove Athena seqs not LANL
  toremove <- c(
    #grep('MK061075',rownames(seqs.all)),
    #          grep('MT570536',rownames(seqs.all)),
    #          grep('EU529856',rownames(seqs.all)),
    #          grep('GQ462061',rownames(seqs.all)),
    #          grep('KU319546',rownames(seqs.all)),
              grep('T901246',rownames(seqs.all)),
              grep('T916739',rownames(seqs.all)),
              grep('T902905',rownames(seqs.all)),
              grep('T919113',rownames(seqs.all)),
              grep('T916159',rownames(seqs.all)),
              grep('T916768',rownames(seqs.all)),
              grep('T911793',rownames(seqs.all)),
              grep('T913982',rownames(seqs.all)),
              grep('T920576',rownames(seqs.all)),
              grep('T918604',rownames(seqs.all)),
              grep('T918200',rownames(seqs.all)),
              grep('T913436',rownames(seqs.all)),
              grep('T920296',rownames(seqs.all))
  )
  toremove <- rownames(seqs.all)[toremove]
  seqs.all <- seqs.all[labels(seqs.all) %notin% toremove,]

# update HXB2 label
#hxb2lab <- 'B.FR.1983.HXB2-LAI-IIIB-BRU.K03455'
hxb2lab <- 'HXB2'
ds[TAXA_L==hxb2lab,TAXA_L:=hxb2lab]
ds[TAXA_L==hxb2lab,SUBTYPE_L:='B']

#db.l			<- unique(subset(ds, SUBTYPE_L %in% c('B','02_AG','C','A1','01_AE','G','D','06_cpx','F1','F'), TAXA_L))
#db.a			<- unique(subset(ds, SUBTYPE %in% c('B','02_AG','C','A1','01_AE','G','D','06_cpx','F1','F'), FASTASampleCode))
db.l			<- unique(subset(ds, select=TAXA_L))
db.a			<- unique(subset(ds, select=FASTASampleCode))
hxb2.lab <- labels(seqs.all)[grepl('K03455',labels(seqs.all))]
keep <- seqs.all[labels(seqs.all) %in% c(db.l[, TAXA_L],db.a[,FASTASampleCode],hxb2.lab),]

do <- data.table(lab=labels(keep))
do[, gp:= 'LANL']
do[grep('Amst_',lab), gp:= 'Ams']
do[grep('NL_',lab), gp:= 'NL']
do <- merge(do,unique(subset(ds,select=c('TAXA_L','SUBTYPE_L'))),by.x='lab',by.y='TAXA_L',all.x=T)

do[!(SUBTYPE_L %in% c('B','02_AG','C','A1','01_AE','G','D','06_cpx','F1','F')),SUBTYPE_L:='Other']
#do[!(SUBTYPE %in% c('B','02_AG','C','A1','01_AE','G','D','06_cpx','F1','F')),SUBTYPE:='Other']

lanl_all <- unique(subset(ds,select=c('TAXA_L','SUBTYPE_L')))
nrow(lanl_all)

# LANL subtypes
lanl <- subset(do,gp=='LANL')
lanl[, Alpha.2.code:= sapply(strsplit(lab,'.',fixed=TRUE),'[[',2),]
lanl <- merge(lanl,dgeo,by='Alpha.2.code',all.x=T)
lanl[is.na(WRLD), WRLD:= 'Unknown']
#lanl[is.na(SUBTYPE_L), SUBTYPE_L:= 'Unknown']
nrow(lanl)
table(lanl$SUBTYPE_L)
(table(lanl$SUBTYPE_L) / sum(table(lanl$SUBTYPE_L)))*100
table(lanl$WRLD,lanl$SUBTYPE_L)
lanl[, WRLD:= factor(WRLD,levels=c('WEurope','CEurope','EEuropeCentralAsia','Asia','Africa','MENA','NorthAm','LaAmCar','Oceania','Unknown'))]
lanl[, SUBTYPE_L:= factor(lanl$SUBTYPE_L,levels=c('B','C','02_AG','01_AE','A1','G','D','F1','06_cpx','Other'))]
table(lanl$WRLD,lanl$SUBTYPE_L)

# get seqs from 'other' non-major subtypes
#db.l.oth			<- unique(subset(ds, SUBTYPE_L %notin% c('B','02_AG','C','A1','01_AE','G','D','06_cpx','F1','F'), TAXA_L))
#db.a.oth			<- unique(subset(ds, SUBTYPE %notin% c('B','02_AG','C','A1','01_AE','G','D','06_cpx','F1','F'), FASTASampleCode))
oth <- seqs.all[labels(seqs.all) %notin% c(db.l[, TAXA_L],db.a[,FASTASampleCode]),]
doth <- data.table(lab=labels(oth))
doth[, gp:= 'LANL']
doth[grep('Amst_',lab), gp:= 'Ams']
doth[grep('NL_',lab), gp:= 'NL']
doth <- merge(doth,unique(subset(ds,select=c('TAXA_L','SUBTYPE_L'))),by.x='lab',by.y='TAXA_L',all.x=T)
lanl.oth <- subset(doth,gp=='LANL')
table(lanl.oth$SUBTYPE_L)
lanl.oth[, Alpha.2.code:= sapply(strsplit(lab,'.',fixed=TRUE),'[[',2),]
lanl.oth <- merge(lanl.oth,dgeo,by='Alpha.2.code',all.x=T)
lanl.oth[is.na(WRLD), WRLD:= 'Unknown']
table(lanl.oth$WRLD)

# athena other subtypes
## add risk groups
doth[, ID:= as.integer(gsub('([A-Za-z]+)_([0-9]+)_([0-9-]+)_([T0-9]+)','\\2',lab))]
doth <- merge(doth,subset(dbas,select=c('PATIENT','MODE')),by.x='ID',by.y='PATIENT',all.x=T)
## add subtypes
doth <- merge(doth,unique(subset(ds,select=c('FASTASampleCode','SUBTYPE'))),by.x='lab',by.y='FASTASampleCode',all.x=T)
ams.oth <- subset(doth,gp=='Ams' & MODE==1)
table(ams.oth$SUBTYPE)

ath.oth <- subset(doth,(gp=='Ams' & MODE!=1) | gp=='NL')
table(ath.oth$SUBTYPE)

# ATHENA subtypes
ath <- subset(do,gp!='LANL')
ath[, ID:= as.integer(gsub('([A-Za-z]+)_([0-9]+)_([0-9-]+)_([T0-9]+)','\\2',lab))]
# add risk groups
ath <- merge(ath,subset(dbas,select=c('PATIENT','MODE')),by.x='ID',by.y='PATIENT',all.x=T)

# add subtypes
ath <- merge(ath,unique(subset(ds,select=c('FASTASampleCode','SUBTYPE'))),by.x='lab',by.y='FASTASampleCode',all.x=T)
ath[!SUBTYPE %in% c('B','02_AG','C','A1','01_AE','G','D','06_cpx','F1','F'), SUBTYPE:='Other']

# Amsterdam MSM
ams <- subset(ath,gp=='Ams' & MODE==1)
table(ams$SUBTYPE)

# Rest of ATHENA
ath_nonAms <- subset(ath,!(gp=='Ams' & MODE==1))
ath_nonAms[!SUBTYPE %in% c('B','02_AG','C','A1','01_AE','G','D','06_cpx','F1','F'), SUBTYPE:='Other']
table(ath_nonAms$SUBTYPE)

# number of Ams non-MSM major subtypes
nrow(subset(ath,(gp=='Ams' & MODE!=1) & SUBTYPE!='Other')) # 1321
# number of non-Ams major subtypes
nrow(subset(ath,gp!='Ams' & SUBTYPE!='Other')) # 7119
# sum = 8440
nrow(subset(ath,!(gp=='Ams' & MODE==1) & SUBTYPE!='Other')) # 8440

# Non-Amsterdam
nl <- subset(ath,gp=='NL')
