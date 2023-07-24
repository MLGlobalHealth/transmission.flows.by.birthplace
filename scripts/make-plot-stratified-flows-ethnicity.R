

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

# stratify by bplace of recipient
po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = 'tpair_prob_w'
)
po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po[, PAIR_ID := as.integer(gsub(paste0('tpair_prob_w\\[([0-9]+)\\]'),'\\1',as.character(variable)))]
tmp <- subset(do, select = c('PAIR_ID','FROM_BPLACE','TO_BPLACE','TRANS_STAGE'))
setnames(tmp,'FROM_BPLACE','BPLACE')
po <- merge(po, tmp, by = 'PAIR_ID')
po <- po[, list(value = sum(value)), by = c('draw','BPLACE','TO_BPLACE')]
tmp <- po[, list(total = sum(value)), by = c('draw','TO_BPLACE')]
po <- merge(po, tmp, by = c('draw','TO_BPLACE'))
po[, paf := value/total]

po <- po[,
         list( q = quantile(paf, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('TO_BPLACE','BPLACE')
]
po <- dcast.data.table(po, TO_BPLACE+BPLACE~stat, value.var = 'q')
setnames(po,'BPLACE','FROM_BPLACE')
po[, L:= paste0(round(M*100,0),paste0('% ['),round(CL*100,0),'-',round(CU*100,0),'%]')]
saveRDS(po,file=paste0(outfile.base,'-rep_',replicate,'-PAF_bplace_stratify','.RDS'))
write.csv(po,file=paste0(outfile.base,'-flows_bplace_stratify','.csv'))
tmp <- subset(po,select=c('TO_BPLACE','FROM_BPLACE','L'))
tmp <- dcast(tmp,TO_BPLACE~FROM_BPLACE,value.var='L')
write.csv(tmp,file=paste0(outfile.base,'-flows_bplace_stratify_table','.csv'))

g <- ggplot(subset(po,TO_BPLACE!='Overall')) + geom_bar(aes(x=TO_BPLACE,y=M,fill=FROM_BPLACE),stat='identity',position=position_dodge(width=0.9)) +
  geom_errorbar(aes(x=TO_BPLACE,ymin=CL, ymax=CU,fill=FROM_BPLACE),position=position_dodge(width=0.9), width=0.5, colour="black")	+
  scale_fill_npg() +
  labs(fill='Birthplace of\nlikely transmitter', y='Contribution to transmissions',x='Birthplace of\nrecipient') +
  theme_bw() +
  theme(legend.pos='bottom') + #,
        #axis.text.x = element_text(angle=0, vjust = 0.5)) + #,
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))

ggsave(file = ('/Users/alexb/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Roadmap/sources/ethnicity_analysis/paper2_figures/flows_frombplace_tobplace_contributions.pdf'), g, w = 8, h = 5)

#g <- ggarrange(g1 + rremove("xlab"),g2+ theme(axis.text.x = element_text(angle=90, vjust = 0.5))+ rremove("xlab") + rremove("ylab"),
#               ncol=2,widths=c(0.35,0.65),align='hv',common.legend=T,legend = "none")
#g_bplace <- annotate_figure(g,bottom = text_grob("Birthplace of recipient",size=28))

# load sampling-adjusted flows
po <- readRDS(file=paste0(outfile.base,'-stratified_flows_frombplace_adjusted_samplingbias','.RDS'))
po <- readRDS(file=paste0(outfile.base,'-adjusted_flows_samplingofcases_bplacecase_bplacesrc','.RDS'))

# plot sankey ----
tmp <- do[, list(N=length(unique(TO_SEQUENCE_ID))),by=c('TO_BPLACE')]
po <- merge(po,tmp,by='TO_BPLACE')
po[, TO_BPLACE:= paste0(gsub('\n',' ',TO_BPLACE),' (N=',N,')')]
po[, FROM_BPLACE:= gsub('\n',' ',FROM_BPLACE)]

nodes <- data.frame(
  name=c(as.character(po$FROM_BPLACE),
         as.character(po$TO_BPLACE)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
po$IDsource <- match(po$FROM_BPLACE, nodes$name)-1
po$IDtarget <- match(po$TO_BPLACE, nodes$name)-1

# Make the Network
pal <- pal_npg('nrc')(7)
my_pal <- paste(shQuote(c(pal,pal), type="cmd"), collapse=", ")
groups <- paste(shQuote(unique(c(as.character(po$FROM_BPLACE),po$TO_BPLACE)), type="cmd"), collapse=", ")
my_color <- paste0('d3.scaleOrdinal() .domain([',groups,']) .range([',my_pal,'])')
my_color <- 'd3.scaleOrdinal() .domain(["Netherlands","W.Europe,\nN.America,Oceania", "Suriname &\nDutch Caribbean",  "S. America &\n Caribbean",
                                         "E. & C. Europe","MENA","Other",
                                         "Netherlands (N=236)","W.Europe, N.America,Oceania (N=48)", "Suriname & Dutch Caribbean (N=34)",
                                          "S. America &  Caribbean (N=40)","E. & C. Europe (N=18)","MENA (N=15)","Other (N=18)"
                                        ]) .range(["#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF",
                                        "#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF"])'

my_color <- 'd3.scaleOrdinal() .domain(["A", "B", "C", "D", "E", "F", "G",
                                        "H", "I", "J", "K", "L", "M", "N"]) .range(["#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF",
                                        "#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF"])'
#my_color <- 'd3.scaleOrdinal() .domain(["A", "B", "C", "D"]) .range(["#E64B35FF", "#4DBBD5FF", "#E64B35FF", "#4DBBD5FF"])'

p <- sankeyNetwork(Links = po, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "M", NodeID = "name",
                   #fontSize= 15, nodeWidth = 20, margin = list(left = 100),
                   sinksRight=FALSE, colourScale=my_color, iterations = 0)
p
onRender(
  p,
  '
  function(el, x) {
    d3.selectAll(".node text").attr("text-anchor", "begin").attr("x", 20);
  }
  '
)

require(htmlwidgets)
#saveWidget(p, file="/Users/alexb/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Roadmap/sources/ethnicity_analysis/paper2_figures/flows_bplace_sankey.html")
saveWidget(p, file="/Users/alexb/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Roadmap/sources/ethnicity_analysis/paper2_figures/flows_frombplace_tobplace_sankey.html")

po[, FROM_BPLACE:= factor(FROM_BPLACE,
                          levels=c('Netherlands','W.Europe, N.America,Oceania','Suriname & Dutch Caribbean',
                                   'S. America &  Caribbean','E. & C. Europe','MENA','Other'))]
po[, TO_BPLACE:= factor(TO_BPLACE,
                          levels=c('Netherlands (N=236)','W.Europe, N.America,Oceania (N=48)','Suriname & Dutch Caribbean (N=34)',
                                   'S. America &  Caribbean (N=40)','E. & C. Europe (N=18)','MENA (N=15)','Other (N=18)'))]
# alluvial plot
ggplot(as.data.frame(subset(po,select=c('FROM_BPLACE','TO_BPLACE','M'))),
       aes(y = M, axis1 = FROM_BPLACE, axis2 = TO_BPLACE)) +
  geom_alluvium(width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Birth place of likely transmitter", "Birth place of incident case"), expand = c(.05, .05)) +
  scale_fill_npg() +
  theme_bw() +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank())
  #scale_fill_brewer(type = "qual", palette = "Set1")

# using ggplot
ggplot(po, aes(x = x,
               next_x = next_x,
               node = node,
               next_node = next_node,
               fill = factor(node),
               label = node)) +
  geom_sankey(flow.alpha = 0.5, node.color = 1) +
  geom_sankey_label(size = 3.5, color = 1, fill = "white") +
  scale_fill_viridis_d() +
  theme_sankey(base_size = 16) +
  theme(legend.position = "none")

pl <- ggplot(df, aes(x = x,
                     next_x = next_x,
                     node = node,
                     next_node = next_node,
                     fill = factor(node),
                     label = node))             # This Creates a label for each node

pl <- pl +geom_sankey(flow.alpha = 0.5,          #This Creates the transparency of your node
                      node.color = "black",     # This is your node color
                      show.legend = TRUE)        # This determines if you want your legend to show

pl <- pl + geom_sankey_label(Size = 3,
                             color = "black",
                             fill = "white") # This specifies the Label format for each node


pl <- pl + theme_bw()
pl <- pl + theme(legend.position = 'none')
pl <- pl + theme(axis.title = element_blank(),
                 axis.text.y = element_blank(),
                 axis.ticks = element_blank(),
                 panel.grid = element_blank())

pl <- pl + scale_fill_viridis_d(option = "inferno")
pl <- pl + labs(title = "Creating a Sankey Diagram")
pl <- pl + labs(subtitle = "Using a simplified ficticious data")
pl <- pl + labs(caption ="Opeyemi Omiwale" )
pl <- pl + labs(fill = 'Nodes')
pl
