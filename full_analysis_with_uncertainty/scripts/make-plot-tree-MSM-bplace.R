library(data.table)
library(tidyverse)
library(tidytree)
library(phyloscannerR)
library(ggtree)
library(pals)
library(glue)
library(ggsci)

indir <- "~/Box Sync/Roadmap/analysis_220713/trees/"
outdir <- '/Users/alexb/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Roadmap/sources/ethnicity_analysis'
trsm <- 'MSM'
rdas <- list.files(indir,pattern = "*AmsMSM__workspace.rda")

# load metadata
pairs.dir <- file.path('~/Documents/GitHub/transmission.flows.by.birthplace','out_Amsterdam','agegps_sensanalysis_210216_MSM-2010_2022')#,'agegps_updated_criteria_MSM-2010_2022')
all_pairs <- readRDS(file=file.path(pairs.dir, 'all_pairs.rds'))
pairs <- readRDS(file.path(pairs.dir,paste0(trsm,"_pairs.rds")))

#for(i in 1:length(rdas)){
 #i <- 7 # Bc4
 #i <- 4 # Bc1
 #i <- 9 # D
 #i <- 1 # 01AE
 #i <- 3 # A1
  #i <- 6 # bc3
  #i <- 5 # bc2
  i <- 2 # 02AG
 a.rda <- rdas[i]
  load(file.path(indir,a.rda))

  subtype <- str_match(a.rda, ".*_subtype_([A-Za-z0-9]+)+_wOutgroup_*")[,2]
  print(subtype)
  ptree <- phyloscanner.trees[[1]]

  tree <- ptree$tree

  read.counts <- ptree$read.counts#[1:6487]

  attr(tree, 'BLACKLISTED') <- c(is.na(read.counts), rep(FALSE, tree$Nnode))

  read.counts <- c(read.counts, rep(1, tree$Nnode))
  read.counts[which(is.na(read.counts))] <- 1

  attr(tree, 'READ_COUNT') <- read.counts

  # load metadata
  d1 <- data.table(tip=tree$tip.label)
  regex.tip.label <- '^([A-Za-z]+)_+(T[0-9]+)_([0-9]+)_([a-zA-Z]+)_([a-zA-Z]+)_([a-zA-Z]+)-([a-zA-Z]+)_([a-zA-Z_]+)_([0-9]+)-([0-9-]+)-([0-9-]+)_([0-9]+)$'
  d1[, ID:= as.numeric(gsub(regex.tip.label,'\\3',tip))]
  d1[, TRSM:= gsub(regex.tip.label,'\\1',tip)]
  d1 <- merge(d1,unique(subset(pairs,select=c('TO_SEQUENCE_ID','TO_AGE','TO_AGE_GP','TO_EST_INFECTION_DATE','TO_COUNTRY'))),by.x='ID',by.y='TO_SEQUENCE_ID',all.x=T)
  d1 <- merge(d1,unique(subset(pairs,select=c('FROM_SEQUENCE_ID','FROM_AGE','FROM_AGE_GP','FROM_EST_INFECTION_DATE','FROM_COUNTRY'))),by.x='ID',by.y='FROM_SEQUENCE_ID',all.x=T)
  d1 <- merge(d1,unique(subset(all_pairs,select=c('TO_SEQUENCE_ID','TO_CLUSTER_NUMBER'))),by.x='ID',by.y='TO_SEQUENCE_ID',all.x=T)
  d1 <- merge(d1,unique(subset(all_pairs,select=c('FROM_SEQUENCE_ID','FROM_CLUSTER_NUMBER'))),by.x='ID',by.y='FROM_SEQUENCE_ID',all.x=T)

  # add birthplace
  d1 <- merge(d1,subset(dind,select=c('PATIENT','ORIGIN')),by.x='ID',by.y='PATIENT',all.x=T)

  # keep first obs for sources that have multiple recipients
  setkey(d1,tip)
  d1 <- d1[,.SD[1],by = tip]

  d1[, prs:= 0]
  d1[TRSM=='AmsMSM', prs := 1]
  d1[ID %in% pairs$FROM_SEQUENCE_ID, prs:= 2]
  d1[ID %in% pairs$TO_SEQUENCE_ID, prs:= 3]
  d1[, bplace:= 0]
  d1[TRSM=='AmsMSM' & ORIGIN!='NL', bplace:= 1]
  d1[TRSM=='AmsMSM' & ORIGIN=='NL', bplace:= 2]

  # flag the individuals who are possible sources to the recipients
  tmp <- unique(subset(d1,select=c('TO_CLUSTER_NUMBER'),prs==3))
  tmp[, poss_src:= 1]
  setnames(tmp,'TO_CLUSTER_NUMBER','FROM_CLUSTER_NUMBER')
  d1 <- merge(d1,tmp,by='FROM_CLUSTER_NUMBER',all.x=T)
  d1[is.na(poss_src), poss_src:= 0]
  #d1[poss_src==T, poss_src_cat:=prs]
  #d1[is.na(poss_src_cat), poss_src_cat:=4]

  x <- full_join(as_tibble(tree), d1, by = c("label" = "tip"))
  tree2 <- as.treedata(x)

  #new.branch.colours <- attr(tree, "BRANCH_COLOURS")
  new.branch.colours <- attr(tree, "BRANCH_COLOURS")
  new.branch.colours <- as.character(new.branch.colours)
  new.branch.colours[new.branch.colours!='AmsMSM'] <- 'non-MSM'
  new.branch.colours <- factor(new.branch.colours, levels = c("AmsMSM",
                                                              "non-MSM"),
                               labels=c('Amsterdam MSM subgraphs','Amsterdam non-MSM or\noutside Amsterdam'))

  new.individual <- tree2@data$prs
  new.individual <- factor(new.individual, levels=c(0,1,2,3),
                           labels=c('Non-Amsterdam MSM','Amsterdam MSM',
                                    'Amsterdam MSM - potential source' ,'Amsterdam MSM - infected 2010-2021'))
  #border <- tree2@data$poss_src
  #border <- factor(border, levels=c(0,1),
  #                 labels=c('Not a recipient or possible source','Recipient or possible source'))
  border <- tree2@data$bplace
  border <- factor(border, levels=c(0,1,2),
                           labels=c('Non-Amsterdam','Born in the Netherlands','Born outside the Netherlands'))

    huh2 <- c(pal_npg("nrc")(1),'grey50')
    pal <- pal_npg("nrc")(9)
    huh2 <- c('grey50',pal[4],pal[2],pal[1])
    pal_branch <- c(pal[3],'grey50')

    tips <- tree$tip.label %>% length()
    #tree.display <- ggtree(tree, color='grey50', size = 0.1) + #%<+% d1 +
      tree.display <- ggtree(tree, aes(color=new.branch.colours), size = 0.3) + #%<+% d1 +
        #scale_color_manual(values = pal_branch, labels=c('Amsterdam MSM subgraphs',
        #                                                 'Amsterdam non-MSM or\noutside Amsterdam'),
        #                   na.value = "black", drop = F, name = "") +
        #geom_point2(aes(color=new.individual), size=0.66) +
      geom_tippoint(aes(color=new.individual,shape=border), size=2) + # 0.75
      #geom_tiplab(size=0.2) +
      #geom_tippoint(aes(color=FROM_AGE_GP), size=0.66) +
      geom_treescale(x=-0.01, y = length(tree$tip.label)*0.5,  offset=4) +
        scale_shape_manual(values = c("Non-Amsterdam" = 20,
                                      "Born in the Netherlands" = 15,
                                      "Born outside the Netherlands" = 17),
                           labels = c("Non-Amsterdam" = "Non-Amsterdam",
                                      "Born in the Netherlands" = "Born in the Netherlands",
                                      "Born outside the Netherlands" = "Born outside the Netherlands"),
                           drop = F, name = "") +
        scale_colour_manual(values = c("Amsterdam MSM subgraphs" = pal[3],
                                     "Amsterdam MSM - infected 2010-2021"= pal[1],
                                     "Amsterdam MSM - potential source" = pal[2],
                                     "Amsterdam MSM" = pal[4],
                                     "Non-Amsterdam MSM" = 'grey50',
                                     'NA' = 'grey50'),
                          labels=c("Amsterdam MSM subgraphs" = 'Amsterdam MSM subgraphs',
                                   "Amsterdam MSM - infected 2010-2021" = 'Amsterdam MSM infected in 2010-2021',
                                   "Amsterdam MSM - potential source" = 'Plausible source of an Amsterdam MSM infection in 2010-2021',
                                   "Amsterdam MSM" = 'Other Amsterdam MSM',
                                   "Non-Amsterdam MSM" = 'Other individuals from outside Amsterdam or a different Amsterdam risk group',
                                   'NA' = ''),
                          na.value = "black", drop = F, name = "") +
      #scale_colour_manual(values = huh2, labels=c('Non-Amsterdam MSM',
      #                                                          'Amsterdam MSM',
      #                                                          'Amsterdam MSM - potential source',
      #                                                          'Amsterdam MSM - infected 2010-2021'),
      #                    na.value = "black", drop = F, name = "") +
                          theme(legend.text = element_text(size=11),legend.position='none') +
      guides(col=guide_legend(ncol=1)) +
      ylim(-1, length(tree$tip.label) +1)

    ggsave(file.path(outdir,glue("type_{subtype}_tree_newlabels_{trsm}_coloursubgraphs_tips_bplace.pdf")), width = 12, height = tips*0.08, limitsize = F)

## change borders
    tips <- tree$tip.label %>% length()
    #tree.display <- ggtree(tree, color='grey50', size = 0.1) + #%<+% d1 +
    tree.display <- ggtree(tree, aes(color=new.branch.colours), size = 0.2, layout='fan') + #%<+% d1 +
      geom_tippoint(aes(color=new.individual,shape=border), size=0.9) +
      geom_treescale(x=-0.01, y = length(tree$tip.label)*0.5,  offset=4) +
      scale_shape_manual(values = c("Non-Amsterdam" = 20,
                                    "Born in the Netherlands" = 15,
                                    "Born outside the Netherlands" = 17),
                         labels = c("Non-Amsterdam" = "Non-Amsterdam",
                                    "Born in the Netherlands" = "Born in the Netherlands",
                                    "Born outside the Netherlands" = "Born outside the Netherlands"),
                         drop = F, name = "") +
      scale_colour_manual(values = c("Amsterdam MSM subgraphs" = pal[3],
                                     "Amsterdam MSM - infected 2010-2021"= pal[1],
                                     "Amsterdam MSM - potential source" = pal[2],
                                     "Amsterdam MSM" = pal[4],
                                     "Non-Amsterdam MSM" = 'grey50',
                                     'NA' = 'grey50'),
                          labels=c("Amsterdam MSM subgraphs" = 'Amsterdam MSM subgraphs',
                                   "Amsterdam MSM - infected 2010-2021" = 'Amsterdam MSM infected in 2010-2021',
                                   "Amsterdam MSM - potential source" = 'Plausible source of an Amsterdam MSM infection in 2010-2021',
                                   "Amsterdam MSM" = 'Other Amsterdam MSM',
                                   "Non-Amsterdam MSM" = 'Other individuals from outside Amsterdam or a different Amsterdam risk group',
                                   'NA' = ''),
                          na.value = "black", drop = F, name = "") +
      theme(legend.text = element_text(size=11),legend.position='bottom') +
      guides(col=guide_legend(ncol=1)) +
      ylim(-1, length(tree$tip.label) +1)

    ggsave(file.path(outdir,glue("type_{subtype}_tree_newlabels_{trsm}_coloursubgraphs_tips_borders_bplace_fan.pdf")), width = 10, height = tips*0.02, limitsize = F)


## add other data to the plot
    # Make the original plot
    #p <- ggtree(tree)

    # generate some random values for each tip label in the data
    #d1 <- data.frame(id=tree$tip.label, val=rnorm(30, sd=3))
    d1 <- data.table(tip=tree$tip.label)
    regex.tip.label <- '^([A-Za-z]+)_+(T[0-9]+)_([0-9]+)_([a-zA-Z]+)_([a-zA-Z]+)_([a-zA-Z]+)-([a-zA-Z]+)_([a-zA-Z_]+)_([0-9]+)-([0-9-]+)-([0-9-]+)_([0-9]+)$'
    d1[, ID:= as.numeric(gsub(regex.tip.label,'\\3',tip))]
    d1[, TRSM:= gsub(regex.tip.label,'\\1',tip)]
    d1 <- merge(d1,unique(subset(pairs,select=c('TO_SEQUENCE_ID','TO_AGE','TO_EST_INFECTION_DATE'))),by.x='ID',by.y='TO_SEQUENCE_ID',all.x=T)
    d1 <- merge(d1,unique(subset(pairs,select=c('FROM_SEQUENCE_ID','FROM_AGE','FROM_EST_INFECTION_DATE'))),by.x='ID',by.y='FROM_SEQUENCE_ID',all.x=T)

    # Make a second plot with the original, naming the new plot "dot",
    # using the data you just created, with a point geom.
    p2 <- facet_plot(tree.display, panel="dot", data=d1, geom=geom_point, aes(x=FROM_AGE2), color='red3',size=0.3)
    p2 + theme_tree2()
    ggsave(file.path(outdir,glue("type_{subtype}_tree_newlabels_{trsm}_panel.pdf")), width = 7, height = tips*0.02, limitsize = F)

    # Make some more data with another random value.
    d2 <- data.frame(id=tree$tip.label, value = abs(rnorm(30, mean=100, sd=50)))

    # Now add to that second plot, this time using the new d2 data above,
    # This time showing a bar segment, size 3, colored blue.
    p3 <- facet_plot(p2, panel='bar', data=d2, geom=geom_segment,
                     aes(x=0, xend=value, y=y, yend=y), size=3, color='blue4')

    # Show all three plots with a scale
    p3 + theme_tree2()
