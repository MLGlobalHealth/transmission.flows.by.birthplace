
## preamble ----
require(data.table)
require(ggplot2)
require(ggsci)
require(scales)
require(grid)
require(ggpubr)
require(ggExtra)
require(cowplot)
require(viridis)

if (0)
{
  args_dir <- list(
    source_dir = '~/Documents/GitHub/transmission.flows.by.birthplace',
    indir = '~/Box\ Sync/Roadmap',
    analysis = 'analysis_220713',
    pairs.dir = '/Users/alexb/Documents/GitHub/transmission.flows.by.birthplace/out_Amsterdam/published_data_MSM-2010_2022',
    outdir = '/Users/alexb/Documents/GitHub/transmission.flows.by.birthplace/out_Amsterdam/mm_bgUnif_piGP_221027b-published_data_MSM-2010_2021-860418',
    clock_model = '/Users/alexb/Box Sync/Roadmap/source_attribution/molecular_clock/hierarchical',
    stanModelFile = 'mm_bgUnif_piGP_221027b',
    job_tag = 'published_data_MSM-2010_2022'
  )
}
args <- args_dir

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

## read stanin
cat('\nReading Stan input data...')
infile.stanin <- list.files(args$outdir, pattern=paste0('_stanin.RData$'), recursive=TRUE)[1]
stopifnot(length(infile.stanin)>0)
stopifnot(length(infile.stanin)<=1)
tmp <- load(file.path(args$outdir, infile.stanin))
stopifnot(c('args','stan_data')%in%tmp)
args$analysis = args_dir$analysis
args$indir = args_dir$indir

tmp <- paste0(outfile.base,'-fitted_stan_model.rds')
cat("\n Read fitted dat ", tmp , "\n")
model_fit <- readRDS(file = tmp)

infile.sampling.prob <- file.path(args_dir$source_dir,'data','sampling_prob_byyear_bplace.RDS')

## plot case-adjusted flows by birthplace of case and source ----
cat(" \n --------------------------------  plot adjusted flows by birthplace -------------------------------- \n")
## % flows from each region PER region of birth of cases (i.e. sum to one per regino of birth of recipient)

cat('\nReading sampling probabilities...')
spy <- readRDS(file=infile.sampling.prob)

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
                                   'S. America &\n Caribbean','E. & C. Europe','MENA','Other'),
                          labels=c('Netherlands','W.Europe,\nN.America & Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\nNon-Dutch Caribbean','E. & C. Europe','MENA','Other'))]
po[, TO_BPLACE:= factor(TO_BPLACE,
                        levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                 'S. America &\n Caribbean','E. & C. Europe','MENA','Other'),
                        labels=c('Netherlands','W.Europe,\nN.America & Oceania','Suriname &\nDutch Caribbean',
                                 'S. America &\nNon-Dutch Caribbean','E. & C. Europe','MENA','Other'))]

po[, L:= paste0(round(M*100,1),'% [',round(CL*100,1),'-',round(CU*100,1),'%]')]
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

# make table for paper
tab <- dcast(po,FROM_BPLACE~TO_BPLACE,value.var='L')
saveRDS(tab,file=paste0(outfile.base,'-table_stratified_flows.RDS'))

## get flows from group a to group b out of total flows ----
cat(" \n --------------------------------  plot flows from group a to group b out of total flows -------------------------------- \n")

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
                                   'S. America &\n Caribbean','E. & C. Europe','MENA','Other'),
                          labels=c('Netherlands','W.Europe,\nN.America & Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\nNon-Dutch Caribbean','E. & C. Europe','MENA','Other'))]
po[, TO_BPLACE:= factor(TO_BPLACE,
                        levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                 'S. America &\n Caribbean','E. & C. Europe','MENA','Other'),
                        labels=c('Netherlands','W.Europe,\nN.America & Oceania','Suriname &\nDutch Caribbean',
                                 'S. America &\nNon-Dutch Caribbean','E. & C. Europe','MENA','Other'))]

po[, L:= paste0(round(M*100,1),'% [',round(CL*100,1),'-',round(CU*100,1),'%]')]

saveRDS(po,file=paste0(outfile.base,'-adjusted_flows_atob_samplingofcases_bplacecase_bplacesrc','.RDS'))

tmp <- expand.grid(FROM_BPLACE=unique(po$FROM_BPLACE),
                   TO_BPLACE=unique(po$TO_BPLACE))
po <- merge(po,tmp,all=T)
po[is.na(M), M:=0]

do[, FROM_BPLACE:= factor(FROM_BPLACE,
                          levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\n Caribbean','E. & C. Europe','MENA','Other'),
                          labels=c('Netherlands','W.Europe,\nN.America & Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\nNon-Dutch Caribbean','E. & C. Europe','MENA','Other'))]
do[, TO_BPLACE:= factor(TO_BPLACE,
                        levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                 'S. America &\n Caribbean','E. & C. Europe','MENA','Other'),
                        labels=c('Netherlands','W.Europe,\nN.America & Oceania','Suriname &\nDutch Caribbean',
                                 'S. America &\nNon-Dutch Caribbean','E. & C. Europe','MENA','Other'))]

breaks_to <- do[, list(N_TO=length(unique(TO_SEQUENCE_ID))),by=c('TO_BPLACE')]
breaks_from <- do[, list(N_FROM=length(unique(FROM_SEQUENCE_ID))),by=c('FROM_BPLACE')]
breaks_to <- breaks_to[order(TO_BPLACE),]
breaks_from <- breaks_from[order(FROM_BPLACE),]

breaks_to[, pos_to:= cumsum(N_TO) - N_TO/2 ]
breaks_from[, pos_from:= cumsum(N_FROM) - N_FROM/2]

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
  theme_bw() +
  theme(axis.text.x=element_text(angle=45,vjust=0.7)) +
  scale_x_continuous(breaks = breaks_to$pos_to, labels = breaks_to$TO_BPLACE, expand = c(0, 0.1)) +
scale_y_continuous(breaks = breaks_from$pos_from, labels = breaks_from$FROM_BPLACE, expand = c(0, 0.1))

## make time-shifting sources plot ----
po <- readRDS(file=paste0(outfile.base,'-adjusted_flows_mwmb_by2years_samplingofcases_mcsamples','.RDS'))

po <- po[,
         list( q = quantile(paf, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('YEAR_GP','FROM_MIGRANT')
]
po <- dcast.data.table(po, YEAR_GP + FROM_MIGRANT ~stat, value.var = 'q')

pal <- pal_npg('nrc')(4)[c(1,4)]

g_srcs <- ggplot(subset(po)) +
  geom_errorbar(aes(x=YEAR_GP,ymin=CL, ymax=CU,fill=FROM_MIGRANT),position=position_dodge(width=0.5),width=0.3, colour="black")	+
  geom_point(aes(x=YEAR_GP,y=M,colour=FROM_MIGRANT),size=2,position=position_dodge(width=0.5)) +
  scale_fill_manual(values=pal) +
  scale_colour_manual(values=pal) +
  labs(x='',col='Birthplace of likely\ntransmitter', y='\nEstimated sources of transmissions\namong Amsterdam MSM in\nAmsterdam transmission chains',col='') +
  theme_bw() +
  theme(legend.pos='right',
        #axis.title.x = element_blank(),
        #axis.text.x = element_text(angle=40, vjust = 0.5)) + #,
        axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9),
        strip.background = element_blank()) +#+ #,
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))

 g <- ggarrange(g_flows + theme_bw(base_size=9) + theme(axis.text.x = element_text(angle=50, vjust = 0.95,hjust = 0.9)),
                ggarrange(g_hmap + theme_bw(base_size=9) + theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust = 0.9)),
                          NULL,widths=c(0.95,0.05)),g_srcs + theme_bw(base_size=9),
                ncol=1,nrow=3,align='v',labels='AUTO',font.label=list(size=12),heights=c(0.33,0.42,0.25))

ggsave(file = paste0(outfile.base,'-stratified_flows_heatmap_temporal_labs.pdf'),
       g, w = 7, h = 11)
ggsave(file = paste0(outfile.base,'-stratified_flows_heatmap_temporal_labs.png'),
       g, w = 7, h = 11)
