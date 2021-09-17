# LDSC
# Getting Started  conda安装环境
git clone https://github.com/bulik/ldsc.git
cd ldsc
conda env create --file environment.yml
source activate ldsc

./ldsc.py -h
./munge_sumstats.py -h

# name: ldsc
# channels:
# - bioconda
# dependencies:
# - python=2.7
# - bitarray=0.8
# - nose=1.3
# - pybedtools=0.7
# - pip
# - pip:
#   - scipy==0.18
#   - pandas==0.20
#   - numpy==1.16
# ~
# ~

# -------------------------------------------------------------------------
#part 1  LD Score Estimation Tutorial
# non-partitioned LD Scores
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1kg_eur.tar.bz2
tar -jxvf 1kg_eur.tar.bz2

cd 1kg_eur

python ../ldsc.py --bfile 22 --l2 --ld-wind-cm 1 --out 22

# produce four files:
# 22.log
# 22.l2.M
# 22.l2.M_5_50
# 22.l2.ldscore.gz


# Partitioned LD Scores
# Step 1: Creating an annot file
	python make_annot.py \
		--gene-set-file GTEx_Cortex.GeneSet \
		--gene-coord-file ENSG_coord.txt \
		--windowsize 100000 \
		--bimfile ./1000G_EUR_Phase3_plink/1000G.EUR.QC.22.bim \
		--annot-file GTEx_Cortex.annot.gz

# GTEx_Cortex.GeneSet ENSG_coord.txt 1000G_EUR_Phase3_plink 为下载内容 输出GTEx_Cortex.annot.gz

# or
	# python make_annot.py \
	# 	--bed-file Brain_DPC_H3K27ac.bed \
	# 	--bimfile 1000G.EUR.QC.22.bim \
	# 	--annot-file Brain_DPC_H3K27ac.annot.gz

# Step 2: Computing LD scores with an annot file.
	python ldsc.py \
		--l2 \
		--bfile ./1000G_EUR_Phase3_plink/1000G.EUR.QC.22 \
		--ld-wind-cm 1 \
		--annot GTEx_Cortex.annot.gz \
		--thin-annot \
		--out GTEx_Cortex \
		--print-snps ./hapmap3_snps/hm.22.snp
# Writing LD Scores for 17261 SNPs to GTEx_Cortex.l2.ldscore.gz


#----------------------------------------------------------------------------------------------------------------------------
# part 2  Heritability and Genetic Correlation
# genetic correlation between schizophrenia and bipolar disorder

wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
tar -jxvf eur_w_ld_chr.tar.bz2
# Downloading Summary Statistics
wget www.med.unc.edu/pgc/files/resultfiles/pgc.cross.bip.zip
wget www.med.unc.edu/pgc/files/resultfiles/pgc.cross.scz.zip
#地址失效，在figshare下载

# head pgc.cross.BIP11.2013-05.txt

# snpid hg18chr bp a1 a2 or se pval info ngt CEUaf
# rs3131972	1	742584	A	G	1.092	0.0817	0.2819	0.694	0	0.16055
# rs3131969	1	744045	A	G	1.087	0.0781	0.2855	0.939	0	0.133028
# rs3131967	1	744197	T	C	1.093	0.0835	0.2859	0.869	0	.
# rs1048488	1	750775	T	C	0.9158	0.0817	0.2817	0.694	0	0.836449
# rs12562034	1	758311	A	G	0.9391	0.0807	0.4362	0.977	0	0.0925926
# rs4040617	1	769185	A	G	0.9205	0.0777	0.2864	0.98	0	0.87156
# rs28576697	1	860508	T	C	1.079	0.2305	0.7423	0.123	0	0.74537
# rs1110052	1	863421	T	G	1.088	0.2209	0.702	0.137	0	0.752294
# rs7523549	1	869180	T	C	1.823	0.8756	0.4929	0.13	0	0.0137615

# head pgc.cross.SCZ17.2013-05.txt

# snpid hg18chr bp a1 a2 or se pval info ngt CEUaf
# rs3131972	1	742584	A	G	1	0.0966	0.9991	0.702	0	0.16055
# rs3131969	1	744045	A	G	1	0.0925	0.9974	0.938	0	0.133028
# rs3131967	1	744197	T	C	1.001	0.0991	0.9928	0.866	0	.
# rs1048488	1	750775	T	C	0.9999	0.0966	0.9991	0.702	0	0.836449
# rs12562034	1	758311	A	G	1.025	0.0843	0.7716	0.988	0	0.0925926
# rs4040617	1	769185	A	G	0.9993	0.092	0.994	0.979	0	0.87156
# rs4970383	1	828418	A	C	1.096	0.1664	0.5806	0.439	0	0.201835
# rs4475691	1	836671	T	C	1.059	0.1181	0.6257	1.02	0	0.146789
# rs1806509	1	843817	A	C	0.9462	0.1539	0.7193	0.383	0	0.600917

# Munge Data
tar -jxvf eur_w_ld_chr.tar.bz2
unzip -o pgc.cross.bip.zip
unzip -o pgc.cross.scz.zip
bunzip2 w_hm3.snplist.bz2


../munge_sumstats.py \
--sumstats pgc.scz.full.2012-04.txt \
--N 1252901 \
--out scz \
--chunksize 500000 \
--merge-alleles ../w_hm3.snplist


../munge_sumstats.py \
--sumstats pgc.bip.full.2012-04.txt \
--N 2427220 \
--out bip \
--chunksize 500000 \
--merge-alleles ../w_hm3.snplist


# LD Score Regression
../ldsc.py \
--rg scz.sumstats.gz,bip.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out scz_bip

less scz_bip.log

# compute heritability and the LD Score regression
../ldsc.py \
--h2 scz.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out scz_h2

less scz_h2.log


# -------------------------------------------------------------------------------------------------------------------------------------------------
# part 3  Cell type specific analyses
# Download the data
# Choose which set of gene sets to analyze. Options include Multi_tissue_gene_expr, Multi_tissue_chromatin, GTEx_brain, Cahoy, ImmGen, or Corces_ATAC
cts_name=GTEx_brain
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/LDSC_SEG_ldscores/${cts_name}_1000Gv3_ldscores.tgz
tar -xvzf ${cts_name}_1000Gv3_ldscores.tgz

# baseline model and standard regression weights

#Download the LD scores
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baseline_ldscores.tgz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/weights_hm3_no_hla.tgz
tar -xvzf 1000G_Phase3_baseline_ldscores.tgz
tar -xvzf weights_hm3_no_hla.tgz

# download and format the BMI summary statistics
# Download and format the summary statistics, or use your own.
wget https://data.broadinstitute.org/alkesgroup/UKBB/body_BMIz.sumstats.gz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
bunzip2 w_hm3.snplist.bz2

python munge_sumstats.py \
--out UKBB_BMI \
--merge-alleles w_hm3.snplist \
--chunksize 500000 \
--sumstats body_BMIz.sumstats.gz

# python munge_sumstats.py \
# --sumstats body_BMIz.sumstats.gz \
# --merge-alleles w_hm3.snplist \
# --out UKBB_BMI

# Run the regressions
cts_name=GTEx_brain
python ldsc.py \
    --h2-cts UKBB_BMI.sumstats.gz \
    --ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. \
    --out BMI_${cts_name} \
    --ref-ld-chr-cts $cts_name.ldcts \
    --w-ld-chr weights_hm3_no_hla/weights.
# Performing regression.
# Results printed to BMI_GTEx_brain.cell_type_results.txt
# Name    Coefficient     Coefficient_std_error   Coefficient_P_value
# Neuron  7.93194788527e-09       3.02894244784e-09       0.00441303625204
# Oligodendrocyte 7.32019970874e-10       3.51868270994e-09       0.417599619801
# Astrocyte       -5.76220451287e-09      2.60400594455e-09       0.98654507806
