
## preamble ----
require(data.table)
require(ggplot2)
require(ggsci)
require(scales)
require(grid)
require(ggpubr)
require(ggExtra)
require(cowplot)
require(dplyr)
require(tidyr)
require(lubridate)
require(here)

if (0)
{
  args_dir <- list(
    source_dir = '~/Documents/GitHub/transmission.flows.by.birthplace',
    #indir = '~/Box\ Sync/Roadmap/source_attribution',
    indir = '~/Box\ Sync/Roadmap',
    analysis = 'analysis_220713',
    pairs.dir = '/Users/alexb/Documents/GitHub/transmission.flows.by.birthplace/out_Amsterdam/published_data_MSM-2010_2022',
    outdir = '/Users/alexb/Documents/GitHub/transmission.flows.by.birthplace/out_Amsterdam/mm_bgUnif_piGP_221027b-published_data_MSM-2010_2021-860418',
    stanModelFile = 'mm_bgUnif_piGP_221027b',
    scenario = 15,
    reps = 1,
    rep = 1,
    simulate_data = T,
    job_tag = 'published_data_MSM-2010_2022',
    undiagnosed = '/Users/alexb/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Roadmap/sources/ethnicity_analysis/undiagnosed/undiagnosed_211102-cohort_2010_2015',
    overwrite = 1
  )
  args <- args_dir
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

## read stanin
cat('\nReading Stan input data...')
infile.stanin <- list.files(args$outdir, pattern=paste0('_stanin.RData$'), recursive=TRUE)[1]
stopifnot(length(infile.stanin)>0)
stopifnot(length(infile.stanin)<=1)
tmp <- load(file.path(args$outdir, infile.stanin))
stopifnot(c('args','stan_data')%in%tmp)

cat(outfile.base)

tmp <- paste0(outfile.base,'-fitted_stan_model.rds')
cat("\n Read fitted dat ", tmp , "\n")
model_fit <- readRDS(file = tmp)

### plot prevalence ----
dp <- readRDS(file=paste0(outfile.base,'-prevalence_CIs','.RDS'))

g_prev <- ggplot(subset(dp)) + geom_bar(aes(x=FROM_BPLACE,y=M,fill=FROM_BPLACE),stat='identity',position=position_dodge(width=0.9)) +
  scale_fill_npg() +
  labs(x='Birthplace of\nincident case',fill='Birthplace of likely\ntransmitter', y='Prevalence by place of birth \n(weighted average 2010-2021)') +
  theme_bw(base_size=28) +
  theme(legend.pos='none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9)) +#+ #,
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))
if(args$overwrite){
  ggsave(file = paste0(outfile.base,'-contribution_to_prevalence.pdf'), g_prev, w = 23, h = 10)
  ggsave(file = paste0(outfile.base,'-contribution_to_prevalence.png'), g_prev, w = 16, h = 10)
}

### plot flows by birthplace ----
cat(" \n --------------------------------  plot flows by birthplace -------------------------------- \n")
spy <- readRDS(file=file.path(args$source_dir,'data','sampling_prob_byyear_bplace','.RDS'))

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
po <- merge(po, spy,
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
if(args$overwrite){saveRDS(po,file=paste0(outfile.base,'-adjusted_flows_samplingofcases','.RDS'))}

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
        axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9)) +#+ #,
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))
if(args$overwrite){
  ggsave(file = paste0(outfile.base,'-contribution_to_flows.pdf'), g1, w = 23, h = 10)
  ggsave(file = paste0(outfile.base,'-contribution_to_flows.png'), g1, w = 16, h = 10)
}

### plot flows over prevalence contributions ----
cat(" \n --------------------------------  plot flows over prevalence contributions -------------------------------- \n")

po <- readRDS(file=paste0(outfile.base,'-flows_frombplace_adjusted_samplingbias_mcsamples','.RDS'))

po[, FROM_BPLACE:= factor(FROM_BPLACE,
                          levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\n Caribbean','E. & C. Europe','MENA','Other'),
                          labels=c('Netherlands','W.Europe,\nN.America & Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\nNon-Dutch Caribbean','E. & C. Europe','MENA','Other'))]

po <- merge(po, dp, by='FROM_BPLACE')
po[, contr_prev:= paf/M]
po <- po[,
         list( q = quantile(contr_prev, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('FROM_BPLACE')
]
po <- dcast.data.table(po, FROM_BPLACE~stat, value.var = 'q')
po[, TO_BPLACE:= 'Overall']
if(args$overwrite){saveRDS(po,file=paste0(outfile.base,'-flows_prev_frombplace','.RDS'))}

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
  ggsave(file = paste0(outfile.base,'-contribution_to_flows_over_prev.pdf'), g2, w = 23, h = 10)
  ggsave(file = paste0(outfile.base,'-contribution_to_flows_over_prev.png'), g2, w = 16, h = 10)
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
  ggsave(file = paste0(outfile.base,'-prev_flows_frombplace_contributions_panel.pdf'),
       g, w = 10, h = 5)
  ggsave(file = paste0(outfile.base,'-prev_flows_frombplace_contributions_panel.png'),
       g, w = 10, h = 5)
}
