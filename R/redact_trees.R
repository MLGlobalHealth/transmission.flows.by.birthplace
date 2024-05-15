require(pals)
require(ggtree)
require(glue)

hivc.db.Date2numeric<- function( x )
{
  if(!class(x)%in%c('Date','character'))	return( x )
  x	<- as.POSIXlt(x)
  tmp	<- x$year + 1900
  x	<- tmp + round( x$yday / ifelse((tmp%%4==0 & tmp%%100!=0) | tmp%%400==0,366,365), d=3 )
  x
}

## ---- rmd.chunk.roadmap.200203.redact.trees----
# Ancestral state reconstruction
roadmap.200323.redact.trees <- function(analysis)
{
  require(big.phylo)
  # Need to use edited scripts from big.phylo due to update in treedater command
  #source('/rds/general/project/ratmann_roadmap_data_analysis/live/R/bigphylo functions.R')
  #source("/rds/general/project/ratmann_roadmap_data_analysis/live/R/bigphylo.cmd.2.R")
  require(data.table)
  require(tidyverse)
  require(ape)
  require(phytools)
  require(treedater)
  require(phyloscannerR)

  analysis <- 'analysis_220713'

  home <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
  infile.seqinfo <- file.path(home,'Data','data_220331','SHM_2201_ROADMAP_220331_tblLAB_seq.rda')
  indir.phsc  <- file.path(home,analysis,'phyloscanner')
  infile.georeg <- file.path(home,'misc','NEWGEO_220713.csv')
  outdir <- file.path(home,analysis,'redacted_trees')
  infiles.phsc <- data.table(F=list.files(indir.phsc, pattern='workspace.rda$', full.names=TRUE, recursive=TRUE))
  infiles.phsc[, BASENAME:= gsub('workspace.rda$','',basename(F))]
  infiles.phsc[, REP:= gsub('.*_wOutgroup_([0-9]+)_.*','\\1',basename(F))]
  infiles.phsc[, SELECT:= gsub('.*_rerooted_([A-Za-z]+)_.*','\\1',basename(F))]

  tmp <- data.table(F_TREE=list.files(outdir, pattern='redacted_trees.newick$', full.names=TRUE, recursive=TRUE))
  tmp[, BASENAME:= gsub('redacted_trees.newick$','',basename(F_TREE))]
  infiles.phsc <- merge(infiles.phsc, tmp, by='BASENAME', all.x=TRUE)
  infiles.phsc <- subset(infiles.phsc,is.na(infiles.phsc$F_TREE) & REP=='000' & SELECT=='AmsMSM')
  #
  #   for each tree:
  #   redact personal data
  #
  for(i in seq_len(nrow(infiles.phsc)))
  {
    #   i<- 4
    cat('\nProcess',i)
    infile <- infiles.phsc[i,F]
    trsm <- infiles.phsc[i,SELECT]
    load(infile)

    subtype <- str_match(infile, ".*_subtype_([A-Za-z0-9]+)+_wOutgroup_*")[,2]
    print(subtype)

    ph <- phyloscanner.trees[[1]][['tree']]
    stopifnot( !any( ph$tip.label=='' ) )
    #stopifnot( is.binary(ph) )
    #
    #   drop tips without sequence dates,
    #   while conserving the ancestral state reconstructions
    #

    # get tip label order
    #is_tip <- ph$edge[,2] <= length(ph$tip.label)
    #ordered_tips <- ph$edge[is_tip, 2]
    #ph$tip.label[ordered_tips]

    p <- ggtree(ph) + geom_tiplab()
    labs <- get_taxa_name(p)

    #   extract taxa names
    regex.tip.label <- '^([A-Za-z]+)_+(T[0-9]+)_([0-9]+)_([a-zA-Z\\-]+)_([a-zA-Z]+)_([a-zA-Z]+)-([a-zA-Z]+)_([a-zA-Z_]+)_([0-9]+)-([0-9-]+)-([0-9-]+)_([0-9]+)$'
    #dsubgraphtaxa[, ID:= as.numeric(gsub(regex.tip.label,'\\3',TAXA))]
    #dsubgraphtaxa[, SEQ_YEAR:= as.numeric(gsub(regex.tip.label,'\\9',TAXA))]
    #dsubgraphtaxa[, FULL_NAME:= paste0(ST,'_',REP,'_',SELECT,'_',NAME)]

    dph <- data.table(  TAXA_LABEL=ph$tip.label,
                        TRSM=trsm,
                        LOC=gsub('([A-Za-z]+)___.*','\\1',ph$tip.label),
                        TAXA= gsub('^[^_]+___(.*)','\\1',ph$tip.label),
                        SID_PID = gsub('^[^_]+___([^_]*_[^_]*)_(.*)','\\1',ph$tip.label),
                        PID = gsub('^[^_]+___[^_]*_([^_]*)_(.*)','\\1',ph$tip.label),
                        #YEAR_SAMPLE = gsub('^[^_]+___[^_]*_[^_]*_[A-Za-z_+]_([0-9+])-(.*)','\\1',ph$tip.label),
                        #YEAR_SAMPLE2 = gsub('^[^_]+___([^\\.]+)_([^\\.]+)_([^\\.]+)_([^\\.]+)_([^\\.]+)-([^\\.]+)_([^\\.]+)_([0-9])-(.*)','\\8',ph$tip.label),
                        YEAR_SAMPLE = as.numeric(gsub(regex.tip.label,'\\9',ph$tip.label)),
                        ORIGIN = gsub('^[^_]+___([^_]*_[^_]*)_([A-Za-z-]*)_([A-Z]*)_(.*)','\\3',ph$tip.label),
                        TAXA_ID= seq_along(ph$tip.label))

    # merge in infection dates
    dinf <- data.table(read.csv(file.path(home,'transmission_sources','Infection_date_est_rec.csv')))
    setnames(dinf,c("id",'estsctodiagMedian','estsctodiagLL','estsctodiagUL'),c("TO_SEQUENCE_ID",'SER_TO_DIAG','SER_TO_DIAG_LL','SER_TO_DIAG_UL'))
    dinf <- unique(dinf)
    dph[LOC %in% c('AmsHSX','AmsMSM','NL','AmsnonHSX','AmsnonMSM'), PATIENT:= PID]
    dph[, PATIENT:= as.integer(PATIENT)]

    dph <- merge(dph,dinf,by.x='PATIENT',by.y='TO_SEQUENCE_ID',all.x=T)
    # calculate infection date
    dph[,DIAGNOSIS_DATE:= as.Date(hiv_pos_d,format="%Y-%m-%d")]
    dph[,DIAGNOSIS_DATE_N:= hivc.db.Date2numeric(DIAGNOSIS_DATE)]
    dph[,EST_INF_DATE:= DIAGNOSIS_DATE_N-SER_TO_DIAG]
    dph[,EST_INF_DATE:= format(date_decimal(EST_INF_DATE), "%Y-%m-%d")]
    dph[,EST_INF_DATE_CL:= DIAGNOSIS_DATE_N-SER_TO_DIAG_UL]
    dph[,EST_INF_DATE_CL:= format(date_decimal(EST_INF_DATE_CL), "%Y-%m-%d")]
    dph[,EST_INF_DATE_CU:= DIAGNOSIS_DATE_N-SER_TO_DIAG_LL]
    dph[,EST_INF_DATE_CU:= format(date_decimal(EST_INF_DATE_CU), "%Y-%m-%d")]

    # read in birth years
    data_age <- readRDS(file.path(home,'transmission_sources','data_age.rds'))
    data_age <- as.data.table(unique(data_age))
    setnames(data_age,c("TO_SEQUENCE_ID","BIRTH_DATE","BIRTH_DATE_DEC"))

    dph <- merge(dph,subset(data_age,select=c('TO_SEQUENCE_ID','BIRTH_DATE')),by.x='PATIENT',by.y='TO_SEQUENCE_ID',all.x=T)
    dph[,BIRTH_DATE:= as.Date(ISOdate(BIRTH_DATE, 1, 1))]
    dph[,EST_INF_DATE:= as.Date(EST_INF_DATE,format="%Y-%m-%d")]
    dph[,AGE_INF:= as.numeric(EST_INF_DATE-BIRTH_DATE)/365]
    dph[, AGE_GP_INF:= cut(AGE_INF,breaks=c(15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100),include.lowest=T,right=F,
                              labels=c('[15-20)','[20-25)','[25-30)','[30-35)','[35-40)','[40-45)','[45-50)','[50-55)','[55-60)','[60-65)','[65-70)',
                                       '[70-75)','[75-80)','[80-85)','[85-90)','[90-95)','[95-100)'))]
    dph[,EST_INF_DATE:= as.character(EST_INF_DATE,format="%Y-%m-%d")]
    dph[,EST_INF_DATE_CU:= as.character(EST_INF_DATE_CU,format="%Y-%m-%d")]
    dph[,EST_INF_DATE_CL:= as.character(EST_INF_DATE_CL,format="%Y-%m-%d")]

    #group birthplaces into coarser strata

    geo <- data.table(read.csv(infile.georeg))
    geo[geo$Alpha_2_code %in%c('AM','AZ','BY','GE','MD','RU','UA'),WRLD:='EEurope']
    geo[geo$Alpha_2_code %in%c('AU','NZ'),WRLD:='Oceania']
    setnames(geo,c('CNTRY','WRLD'),c('CNTRY_born','WRLD_born'))
    dph <- merge(dph,geo,by.x='ORIGIN',by.y='Alpha_2_code',all.x=T)

    dph[, FROM_BPLACE:="Other"]
    dph[WRLD_born %in% c("WEurope","NorthAm","Oceania"), FROM_BPLACE:="W.Europe,\nN.America,Oceania"]
    dph[WRLD_born %in% c("EEurope", "CEurope"), FROM_BPLACE:="E. & C. Europe"]
    dph[WRLD_born %in% c("LaAmCar"), FROM_BPLACE:="S. America &\n Caribbean"]
    dph[WRLD_born %in% c("DutchCarSuriname"), FROM_BPLACE:="Suriname &\nDutch Caribbean"]
    dph[WRLD_born %in% c("MENA"), FROM_BPLACE:="MENA"]
    dph[ORIGIN=="NL", FROM_BPLACE:="Netherlands"]
    # remove spaces and character returns
    dph[, FROM_BPLACE:=gsub('\n','',FROM_BPLACE)]
    dph[, FROM_BPLACE:=gsub(' ','_',FROM_BPLACE)]

    dph[, TAXANEW:= TAXA_LABEL]
    dph[LOC %in% c('AmsHSX','AmsMSM','NL','AmsnonHSX','AmsnonMSM'), TAXANEW:= paste0(LOC,'___',PID,'_',YEAR_SAMPLE)]
    dph[LOC %in% c('AmsHSX','AmsMSM','NL','AmsnonHSX','AmsnonMSM') & is.na(AGE_GP_INF), AGE_GP_INF:='U']
    dph[LOC %in% c('AmsHSX','AmsMSM','NL','AmsnonHSX','AmsnonMSM') & (is.na(EST_INF_DATE) | EST_INF_DATE=='<NA>'), EST_INF_DATE:='U']
    dph[LOC %in% c('AmsHSX','AmsMSM','NL','AmsnonHSX','AmsnonMSM') & is.na(EST_INF_DATE_CU), EST_INF_DATE_CU:='U']
    dph[LOC %in% c('AmsHSX','AmsMSM','NL','AmsnonHSX','AmsnonMSM') & is.na(EST_INF_DATE_CL), EST_INF_DATE_CL:='U']
    dph[LOC %in% c('AmsHSX','AmsMSM','NL','AmsnonHSX','AmsnonMSM'), TAXANEW:= paste0(LOC,'___',AGE_GP_INF,'_',
                                                                                     gsub('-','/',EST_INF_DATE),'_[',gsub('-','/',EST_INF_DATE_CL),
                                                                                     '-',gsub('-','/',EST_INF_DATE_CU),']_',gsub('\n','',FROM_BPLACE))]

    # re-order
    dph <- dph[order(TAXA_ID),]

    # replace tip labels
    ph$tip.label <- dph$TAXANEW

    outfile.phsc <- file.path(outdir, gsub('workspace\\.rda','redacted\\.rda',basename(infile)))
    outfile.tree <- file.path(outdir, gsub('workspace\\.rda','redacted\\.newick',basename(infile)))

    save(ph, file=outfile.phsc)
    write.tree(ph, file=outfile.tree)

    # draw tree
    #   write to file
    pdf(file=gsub('newick','pdf',outfile.tree), w=20, h=10+Ntip(ph)/10)
    plot(ph, show.node.label=TRUE, cex=0.3)
    dev.off()


    outfile <- file.path(outdir,gsub('workspace\\.rda','annotated_tree_redacted.pdf',basename(infile)))
    tmp <- vector('list')
    tmp[['tree']] <- ph
    tmp[['tree']][['node.states']] <- tmp[['tree']][['mapped.edge']] <- tmp[['tree']][['maps']] <- NULL
    attr(tmp[['tree']],'map.order') <- NULL
    attr(tmp[['tree']],'class') <- 'phylo'
    tmp[['read.counts']] <- rep(1, Ntip(ph))
    write.annotated.tree(tmp, outfile, format="pdf", pdf.scale.bar.width = 0.01, pdf.w = 60, pdf.hm = 0.2, verbose = FALSE)

    # plot with nice colours ----

    #ptree <- phyloscanner.trees[[1]]

    #tree <- ptree$tree

    tree <- ph

    #read.counts <- ptree$read.counts[1:6487]

    #attr(tree, 'BLACKLISTED') <- c(is.na(read.counts), rep(FALSE, tree$Nnode))

    #read.counts <- c(read.counts, rep(1, tree$Nnode))
    #read.counts[which(is.na(read.counts))] <- 1

    #attr(tree, 'READ_COUNT') <- read.counts

    new.branch.colours <- attr(tree, "BRANCH_COLOURS")
    new.individual <- attr(tree, "INDIVIDUAL")

    new.branch.colours <- factor(new.branch.colours, levels = c("AmsMSM",
                                                                "AmsnonMSM",
                                                                "NL",
                                                                "WEurope",
                                                                "EEuropeCentralAsia",
                                                                "NorthAm",
                                                                "LaAmCar",
                                                                "Africa",
                                                                "MENA",
                                                                "Asia",
                                                                "Oceania",
                                                                NA))
    new.individual <- factor(new.individual, levels =  c("AmsMSM",
                                                         "AmsnonMSM",
                                                         "NL",
                                                         "WEurope",
                                                         "EEuropeCentralAsia",
                                                         "NorthAm",
                                                         "LaAmCar",
                                                         "Africa",
                                                         "MENA",
                                                         "Asia",
                                                         "Oceania",
                                                         NA))

    huh <- glasbey(n=21)
    huh2 <- c(rev(pal_npg("nrc")(8))[c(1:6,8,7)], huh[c(7,16,21,19)], "grey50")
    tips <- tree$tip.label %>% length()

    tree.display <- ggtree(tree, aes(color=new.branch.colours), size = 0.1) +  geom_tiplab(aes(color=new.individual),size=0.95) +
      geom_point2(aes(color=new.individual), size=0.66) +
      geom_treescale(x=-0.01, y = length(tree$tip.label)*0.5,  offset=4) +
      scale_colour_manual(values = huh2, na.value = "black", drop = F, labels = c("Amsterdam - MSM",
                                                                                  "Amsterdam - non-MSM",
                                                                                  "Netherlands",
                                                                                  "Western Europe",
                                                                                  "Eastern Europe & Central Asia",
                                                                                  "North America",
                                                                                  "Latin America & Caribbean",
                                                                                  "Sub-Saharan Africa",
                                                                                  "Middle-East and North Africa",
                                                                                  "South- and Southeast Asia",
                                                                                  "Oceania",
                                                                                  "Unassigned"), name = "Region/risk group") +
      theme(legend.text = element_text(size=6)) +
      guides(col=guide_legend(ncol=2)) +
      ylim(-1, length(tree$tip.label) +1)
    if(subtype=='Bc4'){
      tree.display <- tree.display + geom_treescale(x=0.01, y = length(tree$tip.label)*0.5,  offset=-150)
    }

    ggsave(file.path(outdir,glue("type_{subtype}_tree_newlabels_{trsm}.pdf")), width = 11, height = tips*0.06, limitsize = F)

  }
}

