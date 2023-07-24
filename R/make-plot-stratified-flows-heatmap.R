
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

### merge in patient metadata ----
load(infile.meta)
dind <- data.table(dind)

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


## plot case-adjusted flows by birthplace of case and source ----
cat(" \n --------------------------------  plot adjusted flows by birthplace -------------------------------- \n")
## % flows from each region PER region of birth of cases (i.e. sum to one per regino of birth of recipient)

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
po <- po[, list(value = sum(value/psi)), by = c('draw','FROM_BPLACE','TO_BPLACE')]
tmp <- po[, list(total = sum(value)), by = c('draw','TO_BPLACE')]
po <- merge(po, tmp, by = c('draw','TO_BPLACE'))
po[, paf := value/total]
po <- po[,
         list( q = quantile(paf, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('FROM_BPLACE','TO_BPLACE')
]
po <- dcast.data.table(po, FROM_BPLACE + TO_BPLACE ~stat, value.var = 'q')
po[, FROM_BPLACE:= factor(FROM_BPLACE,
                          levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\n Caribbean','E. & C. Europe','MENA','Other'))]
po[, TO_BPLACE:= factor(TO_BPLACE,
                        levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                 'S. America &\n Caribbean','E. & C. Europe','MENA','Other'))]

saveRDS(po,file=paste0(outfile.base,'-adjusted_flowsINTO_samplingofcases_bplacecase_bplacesrc','.RDS'))


g_flows <- ggplot(subset(po,TO_BPLACE!='Overall')) + geom_bar(aes(x=TO_BPLACE,y=M,fill=FROM_BPLACE),stat='identity',position=position_dodge(width=0.9)) +
  geom_errorbar(aes(x=TO_BPLACE,ymin=CL, ymax=CU,fill=FROM_BPLACE),position=position_dodge(width=0.9), width=0.5, colour="black")	+
  scale_fill_npg() +
  labs(x='Birthplace of incident case',fill='Birthplace of source', y='\nContribution to incident cases in\nMSM born in each geographic region') +
  theme_bw(base_size=28) +
  theme(legend.pos='none',
        axis.title.x = element_blank(),
        #axis.text.x = element_text(angle=40, vjust = 0.5)) + #,
        axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9)) +#+ #,
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))

ggsave(file = paste0(outfile.base,'-adjusted_flowsINTO_samplingofcases_contributions.pdf'),
       g, w = 11, h = 8)
ggsave(file = paste0(outfile.base,'-adjusted_flowsINTO_samplingofcases_contributions.png'),
       g, w = 11, h = 8)


## get flows from group a to group b out of total flows ----
cat(" \n --------------------------------  plot flows from group a to group b out of total flows -------------------------------- \n")

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
po <- merge(po, tmp, by = 'PAIR_ID')
po <- merge(po, subset(spy,select=c('LOC_BIRTH_POS','YEAR_OF_INF_EST','psi')),
            by.x=c('TO_BPLACE','YEAR_OF_INF_EST'),by.y=c('LOC_BIRTH_POS','YEAR_OF_INF_EST'))
po <- po[, list(value = sum(value/psi)), by = c('draw','FROM_BPLACE','TO_BPLACE')]
tmp <- po[, list(total = sum(value)), by = c('draw')]
po <- merge(po, tmp, by = c('draw'))
po[, paf := value/total]
po <- po[,
         list( q = quantile(paf, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('FROM_BPLACE','TO_BPLACE')
]
po <- dcast.data.table(po, FROM_BPLACE + TO_BPLACE ~stat, value.var = 'q')
po[, FROM_BPLACE:= factor(FROM_BPLACE,
                          levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\n Caribbean','E. & C. Europe','MENA','Other'))]
po[, TO_BPLACE:= factor(TO_BPLACE,
                        levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                 'S. America &\n Caribbean','E. & C. Europe','MENA','Other'))]

saveRDS(po,file=paste0(outfile.base,'-adjusted_flows_atob_samplingofcases_bplacecase_bplacesrc','.RDS'))

tmp <- expand.grid(FROM_BPLACE=unique(po$FROM_BPLACE),
                   TO_BPLACE=unique(po$TO_BPLACE))
po <- merge(po,tmp,all=T)
po[is.na(M), M:=0]

#tmp <- do[, list(N_TO=length(unique(TO_SEQUENCE_ID))),by=c('TO_BPLACE')]
#po <- merge(po,tmp,by='TO_BPLACE')
#tmp <- do[, list(N_FROM=length(unique(FROM_SEQUENCE_ID))),by=c('FROM_BPLACE')]
#po <- merge(po,tmp,by='FROM_BPLACE')

breaks_to <- do[, list(N_TO=length(unique(TO_SEQUENCE_ID))),by=c('TO_BPLACE')]
breaks_from <- do[, list(N_FROM=length(unique(FROM_SEQUENCE_ID))),by=c('FROM_BPLACE')]
breaks_to <- breaks_to[order(TO_BPLACE),]
breaks_from <- breaks_from[order(FROM_BPLACE),]

breaks_to[, pos_to:= cumsum(N_TO) - N_TO/2 ]
breaks_from[, pos_from:= cumsum(N_FROM) - N_FROM/2]

#breaks_to <- tmp_to[, list(TO_BPLACE=TO_BPLACE,
#                           pos_to = 0.5 * (cumsum(N_TO) + lag(cumsum(N_TO), default = 0)))]
#breaks_from <- tmp_from[, list(FROM_BPLACE=FROM_BPLACE,
#                               pos_from = 0.5 * (cumsum(N_FROM) + lag(cumsum(N_FROM), default = 0)))]


po <- merge(po,breaks_to,by='TO_BPLACE')
po <- merge(po,breaks_from,by='FROM_BPLACE')

po[, flows_cat:= cut(M,breaks=c(0,0.01,0.02,0.03,0.05,0.10,0.30,0.50,1),
                     label = c('0-1%','1-2%','2-3%','3-5%','5-10%','10-30%','30-50%','>50%'))]

g_hmap <- ggplot(po, aes(x=pos_to, y=pos_from,fill= flows_cat)) +
  geom_tile(aes(
    width = N_TO,
    height = N_FROM),
     color = "black") +
  scale_fill_viridis(discrete=TRUE,option="magma",direction=-1,begin=0.1,na.value = "white") +
  labs(x='Birth place of incident case\n\n',y='\nBirth place of source',fill='Estimated transmission flows') +
  #scale_x_discrete(expand = c(0,0)) +
  #scale_y_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45,vjust=0.7)) +
  #scale_x_continuous(breaks = breaks_from$pos_from, labels = breaks_from$FROM_BPLACE, expand = c(0, 0.1)) +
  #scale_y_continuous(breaks = breaks_to$pos_to, labels = breaks_to$TO_BPLACE, expand = c(0, 0.1))
  scale_x_continuous(breaks = breaks_to$pos_to, labels = breaks_to$TO_BPLACE, expand = c(0, 0.1)) +
scale_y_continuous(breaks = breaks_from$pos_from, labels = breaks_from$FROM_BPLACE, expand = c(0, 0.1))

ggsave(file = paste0(outfile.base,'-adjusted_flows_atob_samplingofcases_bplacecase_bplacesrc_heatmap.pdf'),
       g_hmap, w = 9, h = 7)
ggsave(file = paste0(outfile.base,'-adjusted_flows_atob_samplingofcases_bplacecase_bplacesrc_heatmap.png'),
       g_hmap, w = 9, h = 7)



## make panel plot ----

g <- ggarrange(g_flows + theme_bw(base_size=9) + theme(axis.text.x = element_text(angle=50, vjust = 0.95,hjust = 0.9)),
               ggarrange(g_hmap + theme_bw(base_size=9) + theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust = 0.9)),
                         NULL,widths=c(0.95,0.05)),
               ncol=1,nrow=2,align='v',labels='AUTO',font.label=list(size=14),heights=c(0.47,0.53))

ggsave(file = paste0(outfile.base,'-adjusted_flows_stratified_heatmap.pdf'),
       g, w = 7, h = 8)
ggsave(file = paste0(outfile.base,'-adjusted_flows_stratified_heatmap.png'),
       g, w = 7, h = 8)

## make sankey plot ----

po[, TO_BPLACE_N:= paste0(gsub('\n',' ',TO_BPLACE),' (N=',N_TO,')')]
po[, FROM_BPLACE:= gsub('\n',' ',FROM_BPLACE)]

# Make the Network
pal <- pal_npg('nrc')(7)

po[, FROM_BPLACE:= factor(FROM_BPLACE,
                          levels=c('Netherlands','W.Europe, N.America,Oceania','Suriname & Dutch Caribbean',
                                   'S. America &  Caribbean','E. & C. Europe','MENA','Other'))]
po[, TO_BPLACE_N:= factor(TO_BPLACE_N,
                        levels=c('Netherlands (N=236)','W.Europe, N.America,Oceania (N=48)','Suriname & Dutch Caribbean (N=34)',
                                 'S. America &  Caribbean (N=40)','E. & C. Europe (N=18)','MENA (N=15)','Other (N=18)'))]

# long data
po_long <- to_lodes_form(data.frame(subset(po,select=c('FROM_BPLACE','TO_BPLACE_N','M'))),
                              key = "BPLACE", value = "Group", id = "Cohort",
                              axes = 1:2)
po_long <- data.table(po_long)
po_long[, BPLACE_ID:= factor(BPLACE,levels=c('FROM_BPLACE','TO_BPLACE_N'),labels=c("Birth place of source", "Birth place of incident case"))]

# alluvial plot
g_sankey <- ggplot(po_long,
       aes( x = BPLACE_ID, stratum = Group, alluvium = Cohort,y = M)) +
  geom_flow(aes(fill=Group),width = 1/12) +
  geom_stratum(aes(fill=Group),width = 1/4) +
  geom_text(stat = "stratum", aes(label = Group),size=3) +
  scale_x_discrete(expand = c(0.2, 0.2)) +
  #scale_x_discrete(limits = c("Birth place of likely transmitter", "Birth place of incident case"), expand = c(.05, .05)) +
  #scale_fill_npg() +
  scale_fill_manual(values=   c('Netherlands' = pal[1], 'W.Europe, N.America,Oceania' = pal[2],
                                'Suriname & Dutch Caribbean' = pal[3],'S. America &  Caribbean' = pal[4],
                                'E. & C. Europe' = pal[5], 'MENA' = pal[6], 'Other' = pal[7],
                                'Netherlands (N=236)' = pal[1], 'W.Europe, N.America,Oceania (N=48)' = pal[2],
                                'Suriname & Dutch Caribbean (N=34)' = pal[3],'S. America &  Caribbean (N=40)' = pal[4],
                                'E. & C. Europe (N=18)' = pal[5], 'MENA (N=15)' = pal[6], 'Other (N=18)' = pal[7])) +
  theme_bw() +
  theme(axis.title=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.text.x=element_text(size=12),plot.margin = margin(0, 0, 0, 0, "cm"),
        legend.position='none',panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
ggsave(file = paste0(outfile.base,'-adjusted_flows_sankey.png'),
       g_sankey, w = 7, h = 7)


# matrix of flows
po <- readRDS(file=paste0(outfile.base,'-adjusted_flows_atob_samplingofcases_bplacecase_bplacesrc','.RDS'))

tmp <- expand.grid(FROM_BPLACE=unique(po$FROM_BPLACE),
                   TO_BPLACE=unique(po$TO_BPLACE))
po <- merge(po,tmp,all=T)
po[is.na(M), M:=0]

g_matrix <- ggplot(po, aes(x=TO_BPLACE, y=FROM_BPLACE,fill=M)) +
  geom_tile(color='grey') +
  geom_text(aes(label = paste0(round(M*100, 2),'%')), size=2.5) +
  scale_fill_viridis(discrete=F,option="magma",direction=-1,begin=0.1,na.value = "white",alpha=0.5) +
  labs(x='Birth place of incident case\n\n',y='\nBirth place of source',fill='Estimated transmission flows') +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust = 0.9),legend.pos='none',
        panel.grid.minor = element_blank(),panel.grid.major = element_blank())


## make sankey panel plot ----

g <- ggarrange(g_flows + theme_bw(base_size=9) + theme(axis.text.x = element_text(angle=50, vjust = 0.95,hjust = 0.9)),
               ggarrange(g_sankey + theme_bw(base_size=9) + theme(axis.title=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),
                                                                    axis.text.x=element_text(size=9),plot.margin = margin(0, 0, 0, 0, "cm"),
                                                                    legend.position='none',text = element_text(size = 6),
                                                                  panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()),
                         g_matrix + theme_bw(base_size=9) +
                           theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust = 0.9),legend.pos='none',
                                 panel.grid.minor = element_blank(),panel.grid.major = element_blank()),
                         widths=c(0.6,0.4),align='v',labels=c('B','C'),font.label=list(size=14),heights=c(0.47,0.53)),
               ncol=1,nrow=2,labels=c('A',NA),font.label=list(size=14),heights=c(0.47,0.53))

ggsave(file = paste0(outfile.base,'-adjusted_flows_stratified_sankey_matrix.pdf'),
       g, w = 10, h = 8)
ggsave(file = paste0(outfile.base,'-adjusted_flows_stratified_sankey_matrix.png'),
       g, w = 10, h = 8)
