require(tidyverse)
require(data.table)
require(ggplot2)
require(abind)
require(rstan)
require(knitr)
require(bayesplot)
require(ggsci)
require(ggpubr)
require(cmdstanr)

# setup args ----
args <- list(
  source_dir = '~/Documents/GitHub/transmission.flows.by.birthplace',
  indir = '~/Box\ Sync/Roadmap/source_attribution',
  outdir = 'out_Amsterdam',
  pairs_dir = 'published_data_MSM-2010_2021',
  job_tag = 'published_data_MSM-2010_2021',
  trsm = 'MSM',
  clock_model = '/Users/alexb/Box Sync/Roadmap/source_attribution/molecular_clock/hierarchical',
  stanModelFile = 'mm_bgUnif_piGP_221027b', # 2D HSGP model
  hmc_stepsize = 0.02,
  hmc_num_samples = 15,
  hmc_num_warmup = 10,
  seed = 42,
  chain = 1,
  #reps = 1,
  #bg = 'unif',
  local = 1,
  time_period='2010-2022',
  m1 = 24,
  m2 = 24,
  B = 1.2
)

## command line parsing if any
args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0)
{
  stopifnot(args_line[[1]]=='-source_dir')
  stopifnot(args_line[[3]]=='-stanModelFile')
  stopifnot(args_line[[5]]=='-chain')
  stopifnot(args_line[[7]]=='-indir')
  stopifnot(args_line[[9]]=='-outdir')
  stopifnot(args_line[[11]]=='-jobtag')
  stopifnot(args_line[[13]]=='-trsm')
  stopifnot(args_line[[15]]=='-cmdstan')
  stopifnot(args_line[[17]]=='-pairs_dir')
  stopifnot(args_line[[19]]=='-clock_model')
  stopifnot(args_line[[21]]=='-time_period')
  stopifnot(args_line[[23]]=='-m1')
  stopifnot(args_line[[25]]=='-m2')
  stopifnot(args_line[[27]]=='-B')
  stopifnot(args_line[[29]]=='-local')

  args <- list()
  args[['source_dir']] <- args_line[[2]]
  args[['stanModelFile']] <- args_line[[4]]
  args[['chain']] <- args_line[[6]]
  args[['indir']] <- args_line[[8]]
  args[['outdir']] <- args_line[[10]]
  args[['job_tag']] <- args_line[[12]]
  args[['trsm']] <- args_line[[14]]
  args[['cmdstan']] <- as.integer(args_line[[16]])
  args[['pairs_dir']] <- args_line[[18]]
  args[['clock_model']] <- args_line[[20]]
  args[['time_period']] <- args_line[[22]]
  args[['m1']] <- as.integer(args_line[[24]])
  args[['m2']] <- as.integer(args_line[[26]])
  args[['B']] <- as.integer(args_line[[28]])
  args[['local']] <- as.integer(args_line[[30]])
}
args

if(args$local==1){
  in.dir <- file.path(args$outdir,args$pairs_dir)
}else{
  in.dir <- file.path(args$outdir,args$job_tag)
}

## set other args
args$file_stanModel <- file.path(args$source_dir, 'stan_model_files',paste0(args$stanModelFile,'.stan'))
tmp <- Sys.getenv("PBS_JOBID")
args$job_id <- ifelse(tmp != '', tmp, as.character(abs(round(rnorm(1) * 1e6))) )
#args$job_dir <- args$outdir
#if(args$local==1)
args$job_dir <- file.path(args$outdir,paste0(args$stanModelFile,'-',args$job_tag,'-',args$job_id))
args$savedata <- TRUE
if (grepl("\\[1\\]", args$job_id))  args$savedata <- TRUE

## make job dir
out.dir <- args$job_dir
#if(args$local==1)
dir.create( out.dir )
outfile.base <- file.path(out.dir, paste0(args$stanModelFile,"-",args$job_tag))
if(args$local==0) outfile.base <- file.path(args$job_dir, paste0(basename(args$job_dir)))
#outfile.base <- file.path(out.dir, paste0(args$stanModelFile,"-",args$job_tag,'-',args$job_id))

## save input args
saveRDS( args, file=file.path(args$job_dir, paste0(basename(args$job_dir), '_args.RDS')))

r <- 1

# load pairs data ----

pairs <- readRDS(file.path(in.dir,paste0(args$trsm,"_pairs_with_metadata.rds")))

tmp <- pairs[, list(N=length(unique(FROM_SEQUENCE_ID))),by='TO_SEQUENCE_ID']

cat(paste0("number of recipients (true pairs) = ", length(unique(pairs$TO_SEQUENCE_ID))))
cat(paste0("mean number of sources per recipient = ", mean(tmp$N)))
cat(paste0("median number of sources per recipient = ", median(tmp$N)))
cat(paste0("number of total pairs = ", nrow(pairs)))

# Prepare data
do <- pairs
do[,LOG_TIME_ELAPSED:= log(TIME_ELAPSED)]

# quantiles of time elapsed
quantile(do$TIME_ELAPSED,probs=c(0.5,0.025,0.975))

# for how many pairs is only one direction possible?
dp <- subset(do,select=c('FROM_SEQUENCE_ID','TO_SEQUENCE_ID'))
dp[, FROM_SEQUENCE_ID:= as.numeric(FROM_SEQUENCE_ID)]
dp[, TO_SEQUENCE_ID:= as.numeric(TO_SEQUENCE_ID)]
dp <- dp[, list(pairs= paste0(min(FROM_SEQUENCE_ID,TO_SEQUENCE_ID),'-', max(FROM_SEQUENCE_ID,TO_SEQUENCE_ID))),by=c('FROM_SEQUENCE_ID','TO_SEQUENCE_ID')]
cp <- subset(dp,select=c('pairs'))
cat(paste0('proportion of IDs in which only one direction is possible: ', round(nrow(unique(cp))/nrow(dp)*100,0),'%'))

# define groupings for model ----
do[,FROM_AGE_INT:= round(FROM_AGE)]
do[,TO_AGE_INT:= round(TO_AGE)]

# get unique time elapsed

tmp <- data.table(TIME_ELAPSED = sort(unique(do$TIME_ELAPSED)))
tmp[, IDX_UNIQUE_TE := seq_len(nrow(tmp))]
do <- merge(do, tmp, by = 'TIME_ELAPSED')
do <- do[order(PAIR_ID),]

# get unique pairs ages
tmp <- data.table(unique(cbind(do$FROM_AGE_INT,do$TO_AGE_INT)))
tmp[, IDX_UNIQUE_PAIR := seq_len(nrow(tmp))]
setnames(tmp, "V1", "FROM_AGE_INT")
setnames(tmp, "V2", "TO_AGE_INT")
do <- merge(do, tmp, by = c("FROM_AGE_INT","TO_AGE_INT"))

# re-order by ages
do <- do[order(FROM_AGE_INT,TO_AGE_INT),]
do <- data.table(do)
do[, PAIR_ID := seq(1,nrow(do)) ]


# load signal cone from clock model ----

dps_clock <- readRDS(file = file.path(in.dir,'clock_quantiles.rds'))

do2 <- copy(do)

pal_3 <- pal_npg("nrc")(4)[c(1,3)]
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
  geom_point(data=do2,aes(x=TIME_ELAPSED,y=GEN_DIST,colour=FROM_BPLACE)) +
  theme_bw()+
  scale_colour_npg() +
  labs(x='\n Time elapsed (in years)',y='\n Genetic distance of the pair \n',colour='Birthplace of probable\ntransmitter')+
  theme(legend.position='bottom',
        legend.text = element_text()) +
  guides(colour = guide_legend(title.position = "top")) +
  scale_y_continuous(expand = c(0,0), breaks=seq(0,0.18,0.02),labels=scales::label_percent(accuracy = 1L)) +
  scale_x_continuous(expand = c(0,0), breaks=seq(0,18,2))
p
ggsave(file=file.path(out.dir,paste0('data_in_signal_cone_',args$time_period,'.pdf')), p, w=7, h=5)
ggsave(file=file.path(out.dir,paste0('data_in_signal_cone_',args$time_period,'.png')), p, w=7, h=5)

# load estimated molecular clock ----
cat(" \n ------------- \n Load quantiles from fitted molecular clock model \n ------------ \n")

cm <- readRDS(file.path(args$clock_model,'clock_model_gamma_hier_220315-stan_fit.rds'))
pd <- cm$draws(inc_warmup = FALSE)
po <- list()
tmp <- pd[,,which(grepl('log_alpha1$',dimnames(pd)[[3]]))]
po$log_alpha1 <- as.vector(tmp)
tmp <- pd[,,which(grepl('log_alpha1_pair_sd',dimnames(pd)[[3]]))]
po$log_alpha1_pair_sd <- as.vector(tmp)
tmp <- pd[,,which(grepl('log_phi$',dimnames(pd)[[3]]))]
po$log_phi <- as.vector(tmp)
tmp <- pd[,,which(grepl('log_phi_pair_sd',dimnames(pd)[[3]]))]
po$log_phi_pair_sd <- as.vector(tmp)

log_alpha1 <- median(po$log_alpha1)
log_alpha1_pair_sd <- median(po$log_alpha1_pair_sd)
log_phi <- median(po$log_phi)
log_phi_pair_sd <- median(po$log_phi_pair_sd)

saveRDS(dps_clock,file = paste0(outfile.base,'-clock_quantiles.rds'))

# define stan data ----

stan_data <- list()
stan_data$N <- nrow(do)
stan_data$P <- length(unique(do$TO_SEQUENCE_ID))
stan_data$S <- length(unique(do$FROM_AGE_GP))
stan_data$x <- do$TIME_ELAPSED
stan_data$y <- do$GEN_DIST
stan_data$pt_idx <- factor(do$TO_SEQUENCE_ID,levels=unique(do$TO_SEQUENCE_ID),labels=seq(1,length(unique(do$TO_SEQUENCE_ID)),1))
tmp <- unique(do, by = 'IDX_UNIQUE_TE')[order(IDX_UNIQUE_TE)]
stan_data$Nux = nrow(tmp)
stan_data$ux = tmp$TIME_ELAPSED
age_vec <- seq(1, 90)
stan_data$a <- max(age_vec)

# posterior medians from clock model
stan_data$log_alpha1 <- log_alpha1
stan_data$log_alpha1_pair_sd <- log_alpha1_pair_sd
stan_data$log_phi <- log_phi
stan_data$log_phi_pair_sd <- log_phi_pair_sd

my.matrix <- matrix(0, nrow=length(unique(do$TO_SEQUENCE_ID)), ncol=nrow(do))
for(i in 1:length(unique(do$TO_SEQUENCE_ID))) {my.matrix[i,stan_data$pt_idx==i]=1}
stan_data$pt_map <- my.matrix

stan_data$M <- c(args$m1,args$m2)
stan_data$M_nD <- stan_data$M[1] * stan_data$M[2]

stan_data$n <- 90
age_idx <- seq.int(1,stan_data$n,1)

age_idx_std <- (age_idx - mean(age_idx))/sd(age_idx)

stan_data$ages <- matrix(data=NA,stan_data$n,2)
stan_data$ages[,1] <- age_idx_std
stan_data$ages[,2] <- age_idx_std

# find all coordinates
ages1_grid = data.table(x = seq.int(1,A,1))
ages1_grid[, x_index := 1:nrow(ages1_grid)]
ages2_grid = data.table(y = seq.int(1,A,1))
ages2_grid[, y_index := 1:nrow(ages2_grid)]

grid = as.data.table( expand.grid(x_index = ages1_grid$x_index,
                                  y_index = ages2_grid$y_index) )
grid = merge(grid, ages1_grid, by = 'x_index')
grid = merge(grid, ages2_grid, by = 'y_index')
cat('the number of entries on the grid, N, is', nrow(grid))

do <- merge(do,grid,by.x=c('FROM_AGE_INT','TO_AGE_INT'),by.y=c('x','y'),all.x=T)

stan_data$coordinates <- matrix(NA,nrow(do),2)
stan_data$coordinates[,1] <- do$x_index
stan_data$coordinates[,2] <- do$y_index

tmp <- data.table(FROM_AGE_STD = sort(unique(age_idx_std)))
tmp[, FROM_AGE_INT := seq_len(nrow(tmp))]
do <- merge(do,tmp,by='FROM_AGE_INT',all.x=T)
do <- do[order(PAIR_ID)]
stan_data$ux_src = tmp$FROM_AGE_STD
stan_data$Nux_src = nrow(tmp)
stan_data$x_2_ux_src = do$FROM_AGE_INT

tmp <- data.table(TO_AGE_STD = sort(unique(age_idx_std)))
tmp[, TO_AGE_INT := seq_len(nrow(tmp))]
do <- merge(do,tmp,by='TO_AGE_INT',all.x=T)
do <- do[order(PAIR_ID)]
stan_data$ux_rec = tmp$TO_AGE_STD
stan_data$Nux_rec = nrow(tmp)
stan_data$x_2_ux_rec = do$TO_AGE_INT

stan_data$L <- rep(args$B*max(do$FROM_AGE_STD,do$FROM_AGE_STD),2)
stan_data$D <- 2

# init values ----
stan_init <- list()
y_mix <- rep(c(0.05,0.1,0.3), 3)
stan_init$y_mix <- y_mix[args$chain]

# save inputs ----
## save image before running Stan
tmp <- names(.GlobalEnv)
tmp <- tmp[!grepl('^.__|^\\.|^model$',tmp)]
save(list=tmp, file=paste0(outfile.base, '_stanin.RData') )

# save stan_data object
#rstan::stan_rdump( names(stan_data), file=paste0(outfile.base, '_cmdstanin.R'), envir=list2env(stan_data))

# save stan.data object
if(args$cmdstan==1){
  #rstan::stan_rdump( names(stan_init), file=file.path(args$job_dir, paste0(basename(args$job_dir), '_cmdstaninit.R')), envir=list2env(stan_init))
  rstan::stan_rdump( names(stan_data), file=file.path(args$job_dir, paste0(basename(args$job_dir), '_cmdstanin.R')), envir=list2env(stan_data))
} else{

  # compile stan model ----
  options(mc.cores = parallel::detectCores())
  sim_mixture_compiled <- cmdstanr::cmdstan_model(args$file_stanModel,
                                                  force_recompile = TRUE,
                                                  include_paths = dirname(args$file_stanModel)
  )


  # run Stan ----
  options(mc.cores=parallel::detectCores())

  # run Stan using cmdstan
  model_fit <- sim_mixture_compiled$sample(
    data = stan_data,
    iter_warmup = 5e2,
    iter_sampling = 2e3,
    refresh = 100,
    parallel_chains = 4,
    chains = 4,
    adapt_delta = 0.9,
    save_warmup = TRUE,
    init = function() list(y_mix=y_mix)
  )

  tmp <- paste0(outfile.base,'-fitted_stan_model.rds')
  cat("\n Save fitted data to file ", tmp , "\n")
  model_fit$save_object(file = tmp)

}
