hivc.db.Date2numeric<- function( x )
{
  if(!class(x)%in%c('Date','character'))	return( x )
  x	<- as.POSIXlt(x)
  tmp	<- x$year + 1900
  x	<- tmp + round( x$yday / ifelse((tmp%%4==0 & tmp%%100!=0) | tmp%%400==0,366,365), d=3 )
  x
}


extract_subgraphs = function(indir.phsc){
  infiles <- data.table(F=list.files(indir.phsc, pattern='rda$', full.names=TRUE, recursive=TRUE))
  infiles <- subset(infiles, grepl('subgraphs_',F))
  #    label by transmission group
  infiles[, SELECT:= gsub('^(.*)subgraphs_([A-Za-z]+).rda','\\2',basename(F))]
  #    ST stores the subtype, and ST_CLADE the large subtree (for ST B only)
  infiles[, ST:= gsub('^.*subtype_([^_]+)_.*\\.rda','\\1',basename(F))]
  infiles[, ST_CLADE:= as.integer(gsub('[^c]*c([0-9]+)','\\1',ST))]
  infiles[, ST:= gsub('([^c]*)c([0-9]+)','\\1',ST)]
  #    ID of bootstrap replicate, 000 for analysis on real alignment
  infiles[, REP:= gsub('^.*wOutgroup_([0-9]+)_.*\\.rda','\\1',basename(F))]
  dsubgraphtaxa <- infiles[, {
    infile <- F
    cat('Process',infile,'\n')
    load(infile)
    if(length(subgraphs)==1){
      subgraph.names <- rep(subgraphs[[1]]$subgraph.name, length(subgraphs[[1]]$tip.label))
      subgraph.taxa <- subgraphs[[1]]$tip.label
      subgraph.parent.state <- subgraphs[[1]]$subgraph.parent.state
    }  else{
      subgraph.names <- unlist(sapply( subgraphs, function(subgraph)  rep(subgraph$subgraph.name, length(subgraph$tip.label)) ))
      subgraph.taxa <- unlist(sapply( subgraphs, function(subgraph)  subgraph$tip.label))
      subgraph.parent.state <- unlist(sapply( subgraphs, function(subgraph)  rep(subgraph$subgraph.parent.state, length(subgraph$tip.label))))
    }
    list(NAME=subgraph.names,
         TAXA= subgraph.taxa,
         ORIGINHOST= subgraph.parent.state
    )
  }, by=c('ST','ST_CLADE','REP','SELECT')]
  #    add meta data from taxa names
  regex.tip.label <- '^([A-Za-z]+)___+(T[0-9]+)_([0-9]+)_([A-Za-z0-9_-]+)'
  dsubgraphtaxa[, ID:= as.numeric(gsub(regex.tip.label,'\\3',TAXA))]
  dsubgraphtaxa[, SEQ_YEAR:= as.numeric(gsub(regex.tip.label,'\\3',TAXA))]
  dsubgraphtaxa[, FULL_NAME:= paste0(ST,'_',REP,'_',SELECT,'_',NAME)]

  dsubgraphtaxa <- unique(dsubgraphtaxa)
  return(dsubgraphtaxa)
}

io_saveRDS <- function(obj, work_dir, out_dir, base_name, check_if_saved_n = 0)
{
  cat('\nSave to file ',file.path(out_dir, base_name),'...')
  tmp <- check_if_saved_n
  repeat
  {
    #   comp_9 = xzfile(tmp, compression = 9)
    # 	saveRDS(fit.gqs, comp_9)
    tryCatch(
      {
        saveRDS(obj, file=file.path(work_dir, base_name))
        if(work_dir!=out_dir)
        {
          file.copy(file.path(work_dir, base_name),
                    file.path(out_dir, base_name),
                    overwrite = TRUE,
                    recursive = FALSE,
                    copy.mode = TRUE,
                    copy.date = TRUE
          )
        }
      }, error = function(err) { warning(err) } )
    if(check_if_saved_n<1)
      break
    check_if_saved <- try(readRDS(file=file.path(out_dir, base_name)))
    if(!'try-error'%in%class(check_if_saved))
      break
    tmp <- tmp-1
    if(tmp<=0)
    {
      stop('Failed to save ',check_if_saved_n,' times')
    }
  }
}
