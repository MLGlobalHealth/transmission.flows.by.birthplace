
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
require(lubridate)

if (1)
{
  args <- list(
    source_dir = '~/Documents/GitHub/source.attr.with.infection.time.fork',
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

# load signal cone from clock model ----

dps_clock <- readRDS(file = file.path(in.dir,'clock_quantiles.rds'))

p.palette <- RColorBrewer::brewer.pal(5,'Oranges')
p.alpha <- 0.7

p <- ggplot(data=dps_clock,aes(x=d_TSeqT)) +
  geom_ribbon(data = dps_clock, aes(ymin = q2.5, ymax = q10), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q90, ymax = q97.5), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q10, ymax = q20), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q80, ymax = q90), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q20, ymax = q30), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q70, ymax = q80), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q30, ymax = q40), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q60, ymax = q70), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q40, ymax = q60), fill = p.palette[5], alpha = p.alpha) +
  geom_point(data=do,aes(x=TIME_ELAPSED,y=GEN_DIST,colour=FROM_BPLACE)) +
  theme_bw(base_size=26)+
  scale_colour_npg() +
  labs(x='\n Time elapsed (in years)',y='\n Genetic distance of the pair \n',colour='Ethnicity of probable transmitter')+
  theme(legend.position='bottom',
        legend.text = element_text(size=18))+
  guides(col = guide_legend(nrow = 2,title.position='top'),by.row=T) +
  scale_y_continuous(expand = c(0,0), breaks=seq(0,0.18,0.02),labels=scales::label_percent(accuracy = 1L)) +
  scale_x_continuous(expand = c(0,0), breaks=seq(0,16,2)) +
  coord_cartesian(xlim=c(0,16))
p
ggsave(file=paste0(outfile.base,'distances_time_elapsed_signalcone.pdf'), p, w=11, h=9)
ggsave(file=paste0(outfile.base,'distances_time_elapsed_signalcone.png'), p, w=11, h=9)

tmp <- do[, list(N=length(unique(TO_SEQUENCE_ID))),by=c('TO_BPLACE')]
po <- merge(po,tmp,by='TO_BPLACE')

p <- ggplot(data=dps_clock,aes(x=d_TSeqT)) +
     geom_ribbon(data = dps_clock, aes(ymin = q2.5, ymax = q10), fill = p.palette[1], alpha = p.alpha) +
     geom_ribbon(data = dps_clock, aes(ymin = q90, ymax = q97.5), fill = p.palette[1], alpha = p.alpha) +
     geom_ribbon(data = dps_clock, aes(ymin = q10, ymax = q20), fill = p.palette[2], alpha = p.alpha) +
     geom_ribbon(data = dps_clock, aes(ymin = q80, ymax = q90), fill = p.palette[2], alpha = p.alpha) +
     geom_ribbon(data = dps_clock, aes(ymin = q20, ymax = q30), fill = p.palette[3], alpha = p.alpha) +
     geom_ribbon(data = dps_clock, aes(ymin = q70, ymax = q80), fill = p.palette[3], alpha = p.alpha) +
     geom_ribbon(data = dps_clock, aes(ymin = q30, ymax = q40), fill = p.palette[4], alpha = p.alpha) +
     geom_ribbon(data = dps_clock, aes(ymin = q60, ymax = q70), fill = p.palette[4], alpha = p.alpha) +
     geom_ribbon(data = dps_clock, aes(ymin = q40, ymax = q60), fill = p.palette[5], alpha = p.alpha) +
     geom_point(data=do,aes(x=TIME_ELAPSED,y=GEN_DIST,colour=FROM_BPLACE)) +
     facet_wrap(~FROM_BPLACE) +
     theme_bw(base_size=26) +
     scale_colour_npg() +
     labs(x='\n Time elapsed (in years)',y='\n Genetic distance of the pair \n',colour='Ethnicity of probable transmitter')+
     theme(legend.position='bottom',
                     legend.text = element_text(size=18),
                     strip.background=element_blank())+
     guides(col = guide_legend(nrow = 2,title.position='top'),by.row=T) +
     scale_y_continuous(expand = c(0,0), breaks=seq(0,0.18,0.02),labels=scales::label_percent(accuracy = 1L)) +
     scale_x_continuous(expand = c(0,0), breaks=seq(0,16,2)) +
     coord_cartesian(xlim=c(0,16))
p
ggsave(file=paste0(outfile.base,'distances_time_elapsed_signalcone_facet.pdf'), p, w=11, h=15)
ggsave(file=paste0(outfile.base,'distances_time_elapsed_signalcone_facet.png'), p, w=11, h=15)

# load flows ----

po <- readRDS(file=paste0(outfile.base,'-rep_1-prob_tpair_w','.RDS'))
set(po,NULL,c('FROM_BPLACE','TO_BPLACE'),NULL)
po <- merge(po,subset(do,select=c('PAIR_ID','FROM_BPLACE','TO_BPLACE')),by='PAIR_ID',all.x=T)
tmp <- do[, list(N=length(TO_SEQUENCE_ID)),by=c('FROM_BPLACE')]
po <- merge(po,tmp,by='FROM_BPLACE')
tmp[, N:= paste0('N = ',N)]

p <- ggplot(data=dps_clock,aes(x=d_TSeqT)) +
  geom_ribbon(data = dps_clock, aes(ymin = q2.5, ymax = q10), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q90, ymax = q97.5), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q10, ymax = q20), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q80, ymax = q90), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q20, ymax = q30), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q70, ymax = q80), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q30, ymax = q40), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q60, ymax = q70), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q40, ymax = q60), fill = p.palette[5], alpha = p.alpha) +
  geom_point(data=po,aes(x=TIME_ELAPSED,y=GEN_DIST,colour=M)) +
  geom_text(
    size    = 5,
    data    = tmp,
    mapping = aes(x = Inf, y = Inf, label = N),
    hjust   = 1.05,
    vjust   = 1.5
  ) +
  facet_wrap(~FROM_BPLACE) +
  theme_bw(base_size=26) +
  scale_colour_gradient(name="Posterior median probabilities\nthat a pair is classified\nby the model as a\ntransmission pair",
                        low = "grey90",
                        high = muted("blue"),
                        breaks = seq(0,1,0.25),
                        n.breaks=4,
                        labels = scales::label_number(accuracy = 0.01),
                        guide=guide_colourbar(direction='horizontal',barwidth=15)) +
  labs(x='\n Time elapsed (in years)',y='\n Patristic distance of the pair \n',colour='Ethnicity of probable transmitter')+
  theme(legend.position='bottom',
        legend.text = element_text(size=18),
        strip.background=element_blank())+
  guides(col = guide_legend(nrow = 2,title.position='top'),by.row=T) +
  scale_y_continuous(expand = c(0,0), breaks=seq(0,0.18,0.02),labels=scales::label_percent(accuracy = 1L)) +
  scale_x_continuous(expand = c(0,0), breaks=seq(0,16,2)) +
  coord_cartesian(xlim=c(0,16)) +
  guides(by.col=T,colour = guide_colorbar(barwidth=15))
p
ggsave(file=paste0(outfile.base,'distances_time_elapsed_signalcone_facet_tpair_samplesize.pdf'), p, w=11, h=15)
ggsave(file=paste0(outfile.base,'distances_time_elapsed_signalcone_facet_tpair_samplesize.png'), p, w=11, h=15)



po <- readRDS(file=paste0(outfile.base,'-tpair_frombplace_stage','.RDS'))
po <- merge(po,subset(do,select=c('PAIR_ID','GEN_DIST')),by='PAIR_ID',all.x=T)
