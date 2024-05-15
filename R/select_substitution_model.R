require(rPlant) # to run fasttree within R
require(big.phylo)
require(data.table)
#source("/rdsgpfs/general/project/ratmann_roadmap_data_analysis/live/bigphylo.cmd.2.R")

analysis <- 'analysis_220713'

home <- '/Users/alexb/Box Sync/Roadmap/Data/data_220331'
home <- '/rds/general/project/ratmann_roadmap_data_analysis/live/analysis_220713'

indir <- file.path(home,'alignments')
indir <- file.path(home,'alignments_bs')
outdir <- file.path(home,'fasttree')

#   get UNIX commands for all trees to be processed
df <- data.table(FIN=list.files(indir, pattern='fasta$',full.names=TRUE))
df <- subset(df, grepl('noDRM',FIN))
df[, ST:= gsub('.*_subtype_([A-Z0-9a-z]+)_.*','\\1',basename(FIN))]
#df[, FOUT:= file.path(outdir,gsub('\\.fasta','_ft.newick',basename(FIN)))]
df[, FOUT:= file.path(outdir,gsub('\\.fasta','_ft_gtr.newick',basename(FIN)))]
#cmds <- vector('list', nrow(df))

for(i in 1:nrow(df)){
  Fasttree(file.name=df[i,'FIN'], file.path="", job.name=NULL, out.name=df[i,'FOUT'], args=NULL, type="DNA", model='GTR', gamma=TRUE, stat=TRUE, print.curl=FALSE,  shared.username=NULL, suppress.Warnings=FALSE)
}


cd /rds/general/project/ratmann_roadmap_data_analysis/live/analysis_220713/select_subst_model
CWD=$(pwd)
echo $CWD
# 01AE
/rdsgpfs/general/user/ablenkin/home/anaconda3/envs/phylo/lib/R/library/big.phylo/ext/FastTree -nt -gtr -gamma -log $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_01AE_wOutgroup_ft_gtr.newick.log $CWD/../alignments_bs/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_01AE_wOutgroup_000.fasta > $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_01AE_wOutgroup_000_ft_gtr.newick
/rdsgpfs/general/user/ablenkin/home/anaconda3/envs/phylo/lib/R/library/big.phylo/ext/FastTree -nt -gamma -log $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_01AE_wOutgroup_ft_JC.newick.log $CWD/../alignments_bs/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_01AE_wOutgroup_000.fasta > $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_01AE_wOutgroup_000_ft_JC.newick

# 02AG
/rdsgpfs/general/user/ablenkin/home/anaconda3/envs/phylo/lib/R/library/big.phylo/ext/FastTree -nt -gtr -gamma -log $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_02AG_wOutgroup_ft_gtr.newick.log $CWD/../alignments_bs/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_02AG_wOutgroup_000.fasta > $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_02AG_wOutgroup_000_ft_gtr.newick
/rdsgpfs/general/user/ablenkin/home/anaconda3/envs/phylo/lib/R/library/big.phylo/ext/FastTree -nt -gamma -log $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_02AG_wOutgroup_ft_JC.newick.log $CWD/../alignments_bs/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_02AG_wOutgroup_000.fasta > $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_02AG_wOutgroup_000_ft_JC.newick

# 06cpx
/rdsgpfs/general/user/ablenkin/home/anaconda3/envs/phylo/lib/R/library/big.phylo/ext/FastTree -nt -gtr -gamma -log $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_06cpx_wOutgroup_ft_gtr.newick.log $CWD/../alignments_bs/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_06cpx_wOutgroup_000.fasta > $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_06cpx_wOutgroup_000_ft_gtr.newick
/rdsgpfs/general/user/ablenkin/home/anaconda3/envs/phylo/lib/R/library/big.phylo/ext/FastTree -nt -gamma -log $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_06cpx_wOutgroup_ft_JC.newick.log $CWD/../alignments_bs/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_06cpx_wOutgroup_000.fasta > $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_06cpx_wOutgroup_000_ft_JC.newick

# A1
/rdsgpfs/general/user/ablenkin/home/anaconda3/envs/phylo/lib/R/library/big.phylo/ext/FastTree -nt -gtr -gamma -log $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_A1_wOutgroup_ft_gtr.newick.log $CWD/../alignments_bs/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_A1_wOutgroup_000.fasta > $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_A1_wOutgroup_000_ft_gtr.newick
/rdsgpfs/general/user/ablenkin/home/anaconda3/envs/phylo/lib/R/library/big.phylo/ext/FastTree -nt -gamma -log $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_A1_wOutgroup_ft_JC.newick.log $CWD/../alignments_bs/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_A1_wOutgroup_000.fasta > $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_A1_wOutgroup_000_ft_JC.newick

# B
/rdsgpfs/general/user/ablenkin/home/anaconda3/envs/phylo/lib/R/library/big.phylo/ext/FastTree -nt -gtr -gamma -log $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_B_wOutgroup_ft_gtr.newick.log $CWD/../alignments_bs/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_B_wOutgroup_000.fasta > $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_B_wOutgroup_000_ft_gtr.newick
/rdsgpfs/general/user/ablenkin/home/anaconda3/envs/phylo/lib/R/library/big.phylo/ext/FastTree -nt -gamma -log $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_B_wOutgroup_ft_JC.newick.log $CWD/../alignments_bs/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_B_wOutgroup_000.fasta > $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_B_wOutgroup_000_ft_JC.newick

# C
/rdsgpfs/general/user/ablenkin/home/anaconda3/envs/phylo/lib/R/library/big.phylo/ext/FastTree -nt -gtr -gamma -log $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_C_wOutgroup_ft_gtr.newick.log $CWD/../alignments_bs/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_C_wOutgroup_000.fasta > $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_C_wOutgroup_000_ft_gtr.newick
/rdsgpfs/general/user/ablenkin/home/anaconda3/envs/phylo/lib/R/library/big.phylo/ext/FastTree -nt -gamma -log $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_C_wOutgroup_ft_JC.newick.log $CWD/../alignments_bs/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_C_wOutgroup_000.fasta > $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_C_wOutgroup_000_ft_JC.newick

# D
/rdsgpfs/general/user/ablenkin/home/anaconda3/envs/phylo/lib/R/library/big.phylo/ext/FastTree -nt -gtr -gamma -log $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_D_wOutgroup_ft_gtr.newick.log $CWD/../alignments_bs/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_D_wOutgroup_000.fasta > $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_D_wOutgroup_000_ft_gtr.newick
/rdsgpfs/general/user/ablenkin/home/anaconda3/envs/phylo/lib/R/library/big.phylo/ext/FastTree -nt -gamma -log $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_D_wOutgroup_ft_JC.newick.log $CWD/../alignments_bs/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_D_wOutgroup_000.fasta > $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_D_wOutgroup_000_ft_JC.newick

# F1
/rdsgpfs/general/user/ablenkin/home/anaconda3/envs/phylo/lib/R/library/big.phylo/ext/FastTree -nt -gtr -gamma -log $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_F1_wOutgroup_ft_gtr.newick.log $CWD/../alignments_bs/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_F1_wOutgroup_000.fasta > $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_F1_wOutgroup_000_ft_gtr.newick
/rdsgpfs/general/user/ablenkin/home/anaconda3/envs/phylo/lib/R/library/big.phylo/ext/FastTree -nt -gamma -log $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_F1_wOutgroup_ft_JC.newick.log $CWD/../alignments_bs/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_F1_wOutgroup_000.fasta > $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_F1_wOutgroup_000_ft_JC.newick

# G
/rdsgpfs/general/user/ablenkin/home/anaconda3/envs/phylo/lib/R/library/big.phylo/ext/FastTree -nt -gtr -gamma -log $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_G_wOutgroup_ft_gtr.newick.log $CWD/../alignments_bs/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_G_wOutgroup_000.fasta > $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_G_wOutgroup_000_ft_gtr.newick
/rdsgpfs/general/user/ablenkin/home/anaconda3/envs/phylo/lib/R/library/big.phylo/ext/FastTree -nt -gamma -log $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_G_wOutgroup_ft_JC.newick.log $CWD/../alignments_bs/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_G_wOutgroup_000.fasta > $CWD/ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_G_wOutgroup_000_ft_JC.newick



# reformat to use with CONSEL
cd '/Users/alexb/OneDrive - Imperial College London/Roadmap/analysis_220713/select_subst_model/'
perl GammaLogToPaup.pl ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_01AE_wOutgroup_ft_gtr.newick.log ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_01AE_wOutgroup_ft_JC.newick.log > comp_JC_GTR.txt

# 01AE
perl /usr/local/bin/GammaLogToPaup.pl ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_01AE_wOutgroup_ft_gtr.newick.log ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_01AE_wOutgroup_ft_JC.newick.log > comp_JC_GTR_01AE.txt
perl /usr/local/bin/GammaLogToPaup.pl ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_01AE_wOutgroup_ft_JC.newick.log ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_01AE_wOutgroup_ft_gtr.newick.log > comp_JC_GTR_01AE.txt
makermt --paup comp_JC_GTR_01AE.txt
consel comp_JC_GTR_01AE
catpv comp_JC_GTR_01AE
catci comp_JC_GTR_01AE

# reading comp_JC_GTR.pv
# rank item    obs     au     np |     bp     pp     kh     sh    wkh    wsh |
#    1    1 -9224.4  1.000  1.000 |  1.000  1.000  1.000  1.000  1.000  1.000 |
#    2    2 9224.4  7e-82  8e-22 |      0      0      0      0      0      0 |

# 02AG
perl /usr/local/bin/GammaLogToPaup.pl ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_02AG_wOutgroup_ft_JC.newick.log ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_02AG_wOutgroup_ft_gtr.newick.log > comp_JC_GTR_02AG.txt
makermt --paup comp_JC_GTR_02AG.txt
consel comp_JC_GTR_02AG
catpv comp_JC_GTR_02AG
catci comp_JC_GTR_02AG

# 06cpx
perl /usr/local/bin/GammaLogToPaup.pl ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_06cpx_wOutgroup_ft_gtr.newick.log ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_06cpx_wOutgroup_ft_JC.newick.log > comp_JC_GTR_06cpx.txt
perl /usr/local/bin/GammaLogToPaup.pl ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_06cpx_wOutgroup_ft_JC.newick.log ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_06cpx_wOutgroup_ft_gtr.newick.log > comp_JC_GTR_06cpx.txt
makermt --paup comp_JC_GTR_06cpx.txt
consel comp_JC_GTR_06cpx
catpv comp_JC_GTR_06cpx
catci comp_JC_GTR_06cpx

# A1
perl /usr/local/bin/GammaLogToPaup.pl ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_A1_wOutgroup_ft_gtr.newick.log ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_A1_wOutgroup_ft_JC.newick.log > comp_JC_GTR_A1.txt
perl /usr/local/bin/GammaLogToPaup.pl ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_A1_wOutgroup_ft_JC.newick.log ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_A1_wOutgroup_ft_gtr.newick.log > comp_JC_GTR_A1.txt
makermt --paup comp_JC_GTR_A1.txt
consel comp_JC_GTR_A1
catpv comp_JC_GTR_A1
catci comp_JC_GTR_A1

# B
perl /usr/local/bin/GammaLogToPaup.pl ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_B_wOutgroup_ft_gtr.newick.log ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_B_wOutgroup_ft_JC.newick.log > comp_JC_GTR_B.txt
perl /usr/local/bin/GammaLogToPaup.pl ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_B_wOutgroup_ft_JC.newick.log ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_B_wOutgroup_ft_gtr.newick.log > comp_JC_GTR_B.txt
makermt --paup comp_JC_GTR_B.txt
consel comp_JC_GTR_B
catpv comp_JC_GTR_B
catci comp_JC_GTR_B

# C
perl /usr/local/bin/GammaLogToPaup.pl ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_C_wOutgroup_ft_gtr.newick.log ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_C_wOutgroup_ft_JC.newick.log > comp_JC_GTR_C.txt
perl /usr/local/bin/GammaLogToPaup.pl ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_C_wOutgroup_ft_JC.newick.log ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_C_wOutgroup_ft_gtr.newick.log > comp_JC_GTR_C.txt
makermt --paup comp_JC_GTR_C.txt
consel comp_JC_GTR_C
catpv comp_JC_GTR_C
catci comp_JC_GTR_C

# D
perl /usr/local/bin/GammaLogToPaup.pl ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_D_wOutgroup_ft_gtr.newick.log ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_D_wOutgroup_ft_JC.newick.log > comp_JC_GTR_D.txt
perl /usr/local/bin/GammaLogToPaup.pl ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_D_wOutgroup_ft_JC.newick.log ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_D_wOutgroup_ft_gtr.newick.log > comp_JC_GTR_D.txt
makermt --paup comp_JC_GTR_D.txt
consel comp_JC_GTR_D
catpv comp_JC_GTR_D
catci comp_JC_GTR_D

# F1
perl /usr/local/bin/GammaLogToPaup.pl ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_F1_wOutgroup_ft_gtr.newick.log ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_F1_wOutgroup_ft_JC.newick.log > comp_JC_GTR_F1.txt
perl /usr/local/bin/GammaLogToPaup.pl ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_F1_wOutgroup_ft_JC.newick.log ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_F1_wOutgroup_ft_gtr.newick.log > comp_JC_GTR_F1.txt
makermt --paup comp_JC_GTR_F1.txt
consel comp_JC_GTR_F1
catpv comp_JC_GTR_F1
catci comp_JC_GTR_F1

# G
perl /usr/local/bin/GammaLogToPaup.pl ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_G_wOutgroup_ft_gtr.newick.log ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_G_wOutgroup_ft_JC.newick.log > comp_JC_GTR_G.txt
perl /usr/local/bin/GammaLogToPaup.pl ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_G_wOutgroup_ft_JC.newick.log ROADMAP_220331_Sequences_LANL_aligned_noDRM_subtype_G_wOutgroup_ft_gtr.newick.log > comp_JC_GTR_G.txt
makermt --paup comp_JC_GTR_G.txt
consel comp_JC_GTR_G
catpv comp_JC_GTR_G
catci comp_JC_GTR_G
