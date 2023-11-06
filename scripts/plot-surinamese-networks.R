
require(data.table)
require(ggplot2)
require(ggsci)
require(ggpubr)
require(igraph)
require(GGally)


## preamble ----
require(data.table)
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


out.dir <- '/Users/alexb/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Roadmap/sources/ethnicity_analysis/paper2_figures'


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
                      #variables = 'den_one'[1]
)
po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po[, PAIR_ID := as.integer(gsub(paste0('tpair_prob_w\\[([0-9]+)\\]'),'\\1',as.character(variable)))]
tmp <- subset(do, select = c('PAIR_ID','FROM_SEQUENCE_ID','TO_SEQUENCE_ID','FROM_BPLACE','TO_BPLACE','TRANS_STAGE'))
setnames(tmp,'FROM_BPLACE','BPLACE')
po <- merge(po, tmp, by = 'PAIR_ID')

po <- po[,
         list( M = quantile(value, probs = c(0.5) )),
         by = c('FROM_SEQUENCE_ID','TO_SEQUENCE_ID')
]
po <- dcast.data.table(po, FROM_SEQUENCE_ID+TO_SEQUENCE_ID~.,value.var='M')
setnames(po,'.','M')

## plot links between everyone in same cluster ----
sg_all <- copy(do)
sg <- subset(sg_all,select=c('FROM_SEQUENCE_ID','TO_SEQUENCE_ID','TO_CLUSTER_NUMBER','FROM_BPLACE','TO_BPLACE'),TO_BPLACE=='Suriname &\nDutch Caribbean')
tmp1 <- subset(sg,select=c('FROM_SEQUENCE_ID','TO_CLUSTER_NUMBER','FROM_BPLACE'))
tmp2 <- subset(sg,select=c('TO_SEQUENCE_ID','TO_CLUSTER_NUMBER','TO_BPLACE'))
setnames(tmp1,c('FROM_SEQUENCE_ID','TO_CLUSTER_NUMBER','FROM_BPLACE'),c('SEQUENCE_ID','CLU','BPLACE'))
setnames(tmp2,c('TO_SEQUENCE_ID','TO_CLUSTER_NUMBER','TO_BPLACE'),c('SEQUENCE_ID','CLU','BPLACE'))
sg <- rbind(tmp1,tmp2)

network <- graph_from_data_frame(d=sg, directed=T)

# change colour of vertices
col = pal_npg("nrc")(5)
V(network)$cluster<- 'grey50'
V(network)$cluster[which(V(network)$name %in% sg$SEQUENCE_ID[sg$BPLACE=='Suriname &\nDutch Caribbean'])]<-col[4]
V(network)$cluster[which(V(network)$name %in% sg$SEQUENCE_ID[sg$BPLACE=='Netherlands'])]<-col[1]
V(network)$color=V(network)$cluster


colrs <- c(col[1],col[2])
l <- layout_nicely(network)

pal <- c('grey50',col[1],col[4])
pdf(file=file.path(out.dir,'networks_surinamese.pdf'),h=30,w=30)
plot(network, layout=layout_nicely, vertex.size=3, edge.color = "grey50",
     edge.arrow.size=1, edge.size=3,
     vertex.label=NA)
legend(x=-1, y=-0.93, c("Other","Born in Suriname",
                        "Born in Netherlands"), pch=21,
       pt.bg=pal, pt.cex=4, cex=3.5, bty="n", ncol=1)
dev.off()

################

## plot linked PAIRS only ----

# load posterior probs of being a pair

sg_all <- copy(do)

sg <- subset(sg_all,select=c('FROM_SEQUENCE_ID','TO_SEQUENCE_ID','TO_CLUSTER_NUMBER','FROM_BPLACE','TO_BPLACE'),TO_BPLACE=='Suriname &\nDutch Caribbean')
sg2 <- merge(sg,po,by=c('FROM_SEQUENCE_ID','TO_SEQUENCE_ID'),all.x=T)

network <- graph_from_data_frame(d=sg2, directed=T)

# change colour of vertices
col = pal_npg("nrc")(5)
V(network)$cluster<- 'grey50'
V(network)$cluster[which(V(network)$name %in% sg$TO_SEQUENCE_ID[sg$TO_BPLACE=='Suriname &\nDutch Caribbean'])]<-col[3]
V(network)$cluster[which(V(network)$name %in% sg$FROM_SEQUENCE_ID[sg$FROM_BPLACE=='Suriname &\nDutch Caribbean'])]<-col[3]
V(network)$cluster[which(V(network)$name %in% sg$FROM_SEQUENCE_ID[sg$FROM_BPLACE=='Netherlands'])]<-col[1]
V(network)$color=V(network)$cluster

#E(network)$width <- E(network)$M*20 # using tpairprob
E(network)$width <- rescale(exp(E(network)$gamma_dens),to=c(0,30)) # using gammdens
#E(network)$newcolor <- E(network)$M


colrs <- c(col[1],col[2])
l <- layout_nicely(network)

pal <- c('grey50',col[3],col[1])
#pdf(file=file.path(out.dir,'networks_surinamese_v2_widths_labs.pdf'),h=30,w=30)
pdf(file=file.path(out.dir,'networks_surinamese_v2_widths_gammadens.pdf'),h=30,w=30)
plot(network, layout=layout_nicely, vertex.size=3,
     #edge.color = E(network)$newcolor,
     edge.color = "grey50",
     edge.arrow.size=1, #edge.size=3,
     vertex.label.cex = 2 ,
     edge.width= E(network)$width,
     vertex.label=NA)
legend(x=-1, y=-0.93, c("Other","Born in Suriname",
                        "Born in Netherlands"), pch=21,
       pt.bg=pal,
       pt.cex=4, cex=3.5, bty="n", ncol=1)
dev.off()


## just plot the cluster with 5 recipients from suriname ----

sg_clu <- subset(sg2,TO_CLUSTER_NUMBER=='B_000_AmsMSM_AmsMSM-SPLIT434_1')
network_clu <- graph_from_data_frame(d=sg_clu, directed=T)

# change colour of vertices
col = pal_npg("nrc")(5)
V(network_clu)$cluster<- 'grey50'
V(network_clu)$cluster[which(V(network_clu)$name %in% sg_clu$TO_SEQUENCE_ID[sg_clu$TO_BPLACE=='Suriname &\nDutch Caribbean'])]<-col[3]
V(network_clu)$cluster[which(V(network_clu)$name %in% sg_clu$FROM_SEQUENCE_ID[sg_clu$FROM_BPLACE=='Suriname &\nDutch Caribbean'])]<-col[3]
V(network_clu)$cluster[which(V(network_clu)$name %in% sg_clu$FROM_SEQUENCE_ID[sg_clu$FROM_BPLACE=='Netherlands'])]<-col[1]
V(network_clu)$color=V(network_clu)$cluster

#E(network_clu)$width <- E(network_clu)$M*20
#E(network_clu)$newcolor <- E(network_clu)$M
E(network_clu)$width <- rescale(exp(E(network_clu)$gamma_dens),to=c(0,30)) # using gammdens


colrs <- c(col[1],col[2])
l <- layout_nicely(network)

pal <- c('grey50',col[3],col[1])
pdf(file=file.path(out.dir,'networks_surinamese_v2_widths_gammadens_clu.pdf'),h=30,w=30)
#png(file=file.path(out.dir,'networks_surinamese_v2_widths.png'),h=30,w=30)
plot(network_clu, layout=layout_nicely, vertex.size=3,
     #edge.color = E(network)$newcolor,
     edge.color = "grey50",
     edge.arrow.size=1, #edge.size=3,
     vertex.label.cex = 2 ,
     edge.width= E(network_clu)$width,
vertex.label=NA)
legend(x=-1, y=-0.93, c("Other","Born in Suriname",
                        "Born in Netherlands"), pch=21,
       pt.bg=pal,
       pt.cex=4, cex=3.5, bty="n", ncol=1)
dev.off()



## sum the weights within people of same ethnicity across the entire set of pairs

all_w <- subset(sg_all,select=c('FROM_SEQUENCE_ID','TO_SEQUENCE_ID','TO_CLUSTER_NUMBER','FROM_BPLACE','TO_BPLACE'))
all_w <- merge(all_w,po,by=c('FROM_SEQUENCE_ID','TO_SEQUENCE_ID'),all.x=T)

tmp <- all_w[, list(sum_w=sum(M[FROM_BPLACE==TO_BPLACE])),by=c('FROM_BPLACE','TO_BPLACE')]
tmp <- subset(tmp,sum_w>0)



## do the same using gamma dens as weights

all_w <- subset(sg_all,select=c('FROM_SEQUENCE_ID','TO_SEQUENCE_ID','TO_CLUSTER_NUMBER','FROM_BPLACE','TO_BPLACE'))
all_w <- merge(all_w,po,by=c('FROM_SEQUENCE_ID','TO_SEQUENCE_ID'),all.x=T)

tmp <- all_w[, list(sum_w=sum(exp(gamma_dens)[FROM_BPLACE==TO_BPLACE])),by=c('FROM_BPLACE','TO_BPLACE')]
tmp <- subset(tmp,sum_w>0)

# need to divide by # individuals otherwise Dutch --> Dutch is always the biggest?

ids_eth1 <- subset(do,select=c('TO_SEQUENCE_ID','TO_BPLACE'))
ids_eth2 <- subset(do,select=c('FROM_SEQUENCE_ID','FROM_BPLACE'))
setnames(ids_eth1,c('TO_SEQUENCE_ID','TO_BPLACE'),c('SEQUENCE_ID','BPLACE'))
setnames(ids_eth2,c('FROM_SEQUENCE_ID','FROM_BPLACE'),c('SEQUENCE_ID','BPLACE'))
ids_eth <- unique(rbind(ids_eth1,ids_eth2))

N_eth <- ids_eth[, list(N=length(unique(SEQUENCE_ID))),by='BPLACE']

tmp <- merge(tmp,N_eth,by.x='FROM_BPLACE',by.y='BPLACE')
tmp[, rel_sum_weights:= sum_w/N]

# Dutch --> Dutch still the highest, followed by Suriname


