## preamble ----
require(data.table)

if (1)
{
  args <- list(
    source_dir = '~/Documents/GitHub/transmission.flows.by.birthplace',
    in_dir = '~/Box\ Sync/Roadmap',
    analysis = 'analysis_220713',
    #out_dir = '~/Documents/GitHub/transmission.flows.by.birthplace/out_Amsterdam/mm_bgUnif_piGP_221027b-agegps_sensanalysis_210216_MSM-618873',
    out_dir = '/Users/alexb/Documents/GitHub/transmission.flows.by.birthplace/out_Amsterdam/mm_bgUnif_piGP_221027b-update_blace_230714_MSM-2010_2021_no_migration_exclusions-322541',
    stanModelFile = 'mm_bgUnif_piGP_221027b',
    undiagnosed = '~/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Roadmap/sources/ethnicity_analysis/undiagnosed/undiagnosed_211102-cohort_2010_2015',
    #job_tag = 'agegps_sensanalysis_210216_MSM-2010_2022'
    job_tag = 'update_blace_230714_MSM-2010_2021_no_migration_exclusions'
  )
}

# set up env variables
cmd2 <- paste0('SCRIPT_DIR=',args$source_dir,'\n',
               'IN_DIR="',args$in_dir,'"\n',
               'OUT_DIR=',args$out_dir,'\n',
               'JOB_TAG=',args$job_tag,'\n',
               'STAN_MODEL_FILE=',args$stanModelFile,'\n',
               'ANALYSIS="',args$analysis,'"\n',
               'OVERWRITE=0\n',
               'UNDIAG_JOB="', args$undiagnosed,'"\n'
)
# check convergence
tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','check-convergence.R'),
              ' -source_dir $SCRIPT_DIR -stanModelFile $STAN_MODEL_FILE -outdir $OUT_DIR -job_tag $JOB_TAG')
cmd2 <- paste0(cmd2,tmp,'\n')
# make plot of cases, sequenced and subtypes
tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','make-panel-plot-nonB-cases.R'),
              ' -source_dir $SCRIPT_DIR -indir $IN_DIR -outdir $OUT_DIR -analysis $ANALYSIS -stanModelFile $STAN_MODEL_FILE -job_tag $JOB_TAG -undiagnosed $UNDIAG_JOB')
cmd2 <- paste0(cmd2,tmp,'\n')
# make plot of contribution of birthplaces to prevalence, viremic MSM and transmission flows
tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','make-plot-contribution-to-prevalence-flows.R'),
              ' -source_dir $SCRIPT_DIR -indir $IN_DIR -outdir $OUT_DIR -analysis $ANALYSIS -undiagnosed $UNDIAG_JOB -overwrite $OVERWRITE')
cmd2 <- paste0(cmd2,tmp,'\n')
# make plot of stratified flows and temporal trends
tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','make-plot-stratified-flows-heatmap-temporal-trends.R'),
              ' -source_dir $SCRIPT_DIR -indir $IN_DIR -outdir $OUT_DIR -analysis $ANALYSIS -undiagnosed $UNDIAG_JOB -overwrite $OVERWRITE')
cmd2 <- paste0(cmd2,tmp,'\n')
# make table of characteristics
tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','make-table-characteristics.R'),
              ' -source_dir $SCRIPT_DIR -indir $IN_DIR -outdir $OUT_DIR -analysis $ANALYSIS -undiagnosed $UNDIAG_JOB')
cmd2 <- paste0(cmd2,tmp,'\n')
# make summarising incident cases and sources
tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','make-table-cases-sequenced-sources.R'),
              ' -source_dir $SCRIPT_DIR -indir $IN_DIR -outdir $OUT_DIR -analysis $ANALYSIS -undiagnosed $UNDIAG_JOB')
cmd2 <- paste0(cmd2,tmp,'\n')

# write submission file
post.processing.file <- file.path(args$out_dir, 'post_processing.sh')
cat(cmd2, file=post.processing.file)
# set permissions
Sys.chmod(post.processing.file, mode='644')


#	run post processing steps
system(post.processing.file)
