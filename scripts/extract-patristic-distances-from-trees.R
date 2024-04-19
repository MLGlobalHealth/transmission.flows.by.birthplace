#install.packages("adephylo")
require(ape)
require(adephylo)
require(ggplot2)
require(data.table)
library(readr)

args <- list(
  in.dir = '/rds/general/project/ratmann_roadmap_data_analysis/live/analysis_220713',
  out.dir = '/rds/general/project/ratmann_roadmap_data_analysis/live/transmission_sources/patristic_distances',
  #out.dir = 'data_Ams/analysis_220713/trees'
  trsm = 'AmsMSM',
  REP = '000'
)

## command line parsing if any
args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0)
{
  stopifnot(args_line[[1]]=='-in.dir')
  stopifnot(args_line[[3]]=='-out.dir')
  stopifnot(args_line[[5]]=='-trsm')
  stopifnot(args_line[[7]]=='-REP')

  args <- list()
  args[['in.dir']] <- args_line[[2]]
  args[['out.dir']] <- args_line[[4]]
  args[['trsm']] <- args_line[[6]]
  args[['REP']] <- args_line[[8]]
}
args

out.dir <- args$out.dir

infiles <- data.table(FIN=list.files(c(file.path(args$in.dir,'phyloscanner'),
                                       file.path(args$in.dir,'phyloscanner_B')), pattern='\\.newick$',full.names=TRUE))
infiles[, ST:= gsub('.*_subtype_([A-Za-z0-9]+)_.*','\\1',basename(FIN))]
infiles[, SELECT:= gsub('^.*_rerooted_([A-Za-z]+)\\.newick$','\\1', basename(FIN))]
infiles[, REP:= gsub('^.*_wOutgroup_([0-9]+)_.*','\\1', basename(FIN))]
infiles <- subset(infiles,SELECT==args$trsm & REP==args$REP) # we only need for one trsm gp

for(s in 1:nrow(infiles)){
  tree <- read.tree(infiles[s,FIN])
  plot(tree,show.node=TRUE)

  # Calculate pairwise distances of tips of the trees #

  distance <- distTips(tree)
  # Transform to a matrix
  distance_matrix <- as.matrix(distance)
  saveRDS(distance,file=file.path(out.dir, paste0('distance_',infiles[s,SELECT],'_',infiles[s,ST],'_',infiles[s,REP],'.rds')))
  saveRDS(distance_matrix,file=file.path(out.dir, paste0('distance_matrix_',infiles[s,SELECT],'_',infiles[s,ST],'_',infiles[s,REP],'.rds')))

  gen_dist_tab <- reshape2::melt(distance_matrix)
  gen_dist_tab <- as.data.table(gen_dist_tab)

  # Remove pairs with 0 distance #
  gen_dist_tab <- gen_dist_tab[-c(which(gen_dist_tab$value==0)),]

  # Add subtype in both tables
  gen_dist_tab[,subtype:=infiles[s,ST]]

  saveRDS(gen_dist_tab,file=file.path(out.dir, paste0('pairwise_dist_',infiles[s,SELECT],'_',infiles[s,ST],'_',infiles[s,REP],'.rds')))
}

gen_dist_tab <- list()
for(s in c('Bc1','Bc2','Bc3','Bc4')){
  gen_dist_tab[[s]] <- readRDS(file=file.path(args$out.dir, paste0('pairwise_dist_',args$trsm,'_',s,'_',args$REP,'.rds')))
}
gen_dist_tabB <- do.call(`rbind`,gen_dist_tab)
gen_dist_tabB <- as.data.table(gen_dist_tabB)
saveRDS(gen_dist_tabB,file=file.path(out.dir, paste0('pairwise_dist_',args$trsm,'_B','_',args$REP,'.rds')))

subtypes <- unique(infiles$ST)
subtypes <- subtypes[!grepl('Bc',subtypes)]
subtypes <- c(subtypes,'B')
gen_dist_tab <- list()
for(s in subtypes){
  gen_dist_tab[[s]] <- readRDS(file=file.path(out.dir, paste0('pairwise_dist_',args$trsm,'_',s,'_',args$REP,'.rds')))
}
# Create a large data table for both trees

gen_dist <- do.call(`rbind`,gen_dist_tab)
gen_dist$Var1 <- as.character(gen_dist$Var1)
gen_dist$Var2 <- as.character(gen_dist$Var2)
# Remove pairs with background sequences
gen_dist <- gen_dist[grep('Amst',gen_dist$Var1),]
gen_dist <- gen_dist[grep('Amst',gen_dist$Var2),]

# Extract the sequences ID from the tree labels
regex.tip.label <- '^([A-Za-z]+)_+([0-9]+)_([0-9-]+)_(T[0-9]+)$'
regex.tip.label2 <- '^([A-Za-z]+)_+(T[0-9]+)_([0-9]+)_([a-zA-Z-]+)_([a-zA-Z0-9-]+)_([a-zA-Z]+)-([a-zA-Z]+)_([a-zA-Z_]+)_([0-9-]+)_([0-9]+)$'
gen_dist[, Var1_id:= as.numeric(gsub(regex.tip.label,'\\2',Var1))]
gen_dist[, Var2_id:= as.numeric(gsub(regex.tip.label,'\\2',Var2))]
gen_dist[, Var1_id2:= as.numeric(gsub(regex.tip.label2,'\\3',Var1))]
gen_dist[, Var2_id2:= as.numeric(gsub(regex.tip.label2,'\\3',Var2))]
gen_dist[is.na(Var1_id), Var1_id:=Var1_id2]
gen_dist[is.na(Var2_id), Var2_id:=Var2_id2]
set(gen_dist,NULL,c('Var1','Var2','Var1_id2','Var2_id2'),NULL)
setnames(gen_dist,c('Var1_id','Var2_id'),c('Var1','Var2'))

setnames(gen_dist,c("Var1","Var2"),c("FROM_SEQUENCE_ID","TO_SEQUENCE_ID"))
setnames(gen_dist,"value","distance")
# Remove duplicate pairs
gen_dist <- unique(gen_dist,by= c("FROM_SEQUENCE_ID","TO_SEQUENCE_ID"))
saveRDS(gen_dist,file=file.path(out.dir, paste0('pairwise_dist_allSTs_',args$trsm,'_',args$REP,'.rds')))

