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
