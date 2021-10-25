library(tidyverse)
library(rhdf5)
library(snow)

file="/share2/pub/zhenggw/zhenggw/brain/l5_all.agg.loom"
h5f <- H5Fopen(file)
exp <- as.data.frame(t(h5f$matrix))
exp$Gene <- h5f$row_attrs$Gene

# Only keep genes with a unique name and tidy data.
exp <- exp %>% add_count(Gene) %>% 
  filter(n==1) %>%
  select(-n) %>%
  gather(key = column,value=Expr,-Gene) %>%
  as_tibble()

# Load gene coordinates and extend upstream and downstream coordinates by 100kb
gene_coordinates <- 
  read_tsv("./NCBI37.3.gene.loc.extendedMHCexcluded",
           col_names = FALSE,col_types = 'cciicc') %>%
  mutate(start=ifelse(X3-100000<0,0,X3-100000),end=X4+100000) %>%
  select(X2,start,end,1) %>% 
  rename(chr="X2", ENTREZ="X1") %>% 
  mutate(chr=paste0("chr",chr))

# Load mouse to human 1to1 orthologs
m2h <- read_tsv("./m2h.txt",col_types = "iccccc") %>% 
  select(musName,entrez) %>%
  rename(Gene=musName) %>% rename(ENTREZ=entrez)


# Add Cell type names
cell_types <- cbind(column=as.character(paste0("V",1:265)),
                    Lvl1=h5f$col_attrs$TaxonomyRank1,
                    Lvl2=h5f$col_attrs$TaxonomyRank2,
                    Lvl3=h5f$col_attrs$TaxonomyRank3,
                    Lvl4=h5f$col_attrs$TaxonomyRank4,
                    Lvl5=h5f$col_attrs$ClusterName,
                    Description=h5f$col_attrs$Description,
                    NCells=h5f$col_attrs$NCells) %>%  
                    as_tibble() %>%
                    mutate(NCells=as.numeric(NCells))

exp_lvl5 <- inner_join(exp,cell_types,by="column") %>% ungroup() %>% rename(Expr_sum_mean=Expr)

# Get average expression at Lvl4
exp_lvl5 <- exp_lvl5 %>% group_by(Lvl4,Gene) %>% summarise(Expr_sum_mean=mean(Expr_sum_mean))

# Write dictonary for cell type names
dic_lvl5 <- select(cell_types,-column,-NCells) %>% select(Lvl1:Lvl4) %>% unique() %>% mutate(Lvl4_LDSC_names=make.names(Lvl4)) 
write_tsv(dic_lvl5,"./dictionary_cell_type_names.Lvl4.txt")


# Remove not expressed genes
not_expressed <- exp_lvl5 %>% 
  group_by(Gene) %>% 
  summarise(total_sum=sum(Expr_sum_mean)) %>% 
  filter(total_sum==0) %>% 
  select(Gene) %>% unique() 

exp_lvl5 <- filter(exp_lvl5,!Gene%in%not_expressed$Gene)


# Each cell type is scaled to the same total number of molecules.
exp_lvl5 <- exp_lvl5 %>% 
  group_by(Lvl4) %>% 
  mutate(Expr_sum_mean=Expr_sum_mean*1e6/sum(Expr_sum_mean))


# specifitiy is defined as the proportion of total expression performed by the cell type of interest (x/sum(x))
exp_lvl5 <- exp_lvl5 %>% 
  group_by(Gene) %>% 
  mutate(specificity=Expr_sum_mean/sum(Expr_sum_mean)) %>% 
  ungroup()


# Only keep genes with 1to1 orthologs
exp_lvl5 <- inner_join(exp_lvl5,m2h,by="Gene")


# Only keep genes that are tested in MAGMA
exp_lvl5 <- inner_join(exp_lvl5,gene_coordinates,by="ENTREZ")


# Get number of genes that represent 10% of the dataset
n_genes <- length(unique(exp_lvl5$ENTREZ))
n_genes_to_keep <- (n_genes * 0.1) %>% round()


# Save expression profile for other processing
save(exp_lvl5,file = "expression.ready.Rdata")

# Get LDSC input top 10%
write_group = function(df,Cell_type) {
  df <- select(df,Lvl4,chr,start,end,ENTREZ)
  dir.create(paste0("LDSC/Bed"), showWarnings = FALSE,recursive = TRUE)
  write_tsv(df[-1],paste0("LDSC/Bed/",make.names(unique(df[1])),".bed"),col_names = F)
return(df)
}

ldsc_bedfile <- function(d,Cell_type){
  d_spe <- d %>% group_by_(Cell_type) %>% top_n(.,n_genes_to_keep,specificity) 
  d_spe %>% do(write_group(.,Cell_type))
}

# Write LDSC input files
# Filter out genes with expression below 1 TPM.
exp_lvl5 %>% filter(Expr_sum_mean>1) %>% ldsc_bedfile("Lvl4")
