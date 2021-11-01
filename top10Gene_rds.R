#get the gene from the .rds dataset for ldsc
#calculate specifity and get the top 10% geneset

library(tidyverse)
library("rhdf5")
library("snow")
library(Seurat)

liver <- readRDS("Liver_cells.rds")
Idents(liver) <- "anno"
liver_avr <- AverageExpression(liver)
liver_expr <- as.data.frame(liver_avr$SCT)
liver_expr <- liver_expr[apply(liver_expr,1,sum)!=0,]

# exp$Gene <- rownames(exp)
# #Only keep genes with a unique name and tidy data.
# exp <- exp %>% add_count(Gene) %>%
#   filter(n==1) %>%
#   select(-n) %>%
#   gather(key = column,value=Expr,-Gene) %>%
#   as.tibble()
top10_function <-function(exp){
exp$Gene <- rownames(exp)
exp <- exp %>% add_count(Gene) %>% 
  filter(n==1) %>%
  select(-n) %>%
  gather(key = column,value=Expr,-Gene) %>%
  as_tibble()

exp <- exp %>%
  group_by(column) %>%
  mutate(Expr_sum_mean=Expr*1e6/sum(Expr))

exp<- exp %>%
  group_by(Gene) %>%
  mutate(specificity=Expr_sum_mean/sum(Expr_sum_mean)) %>%
  ungroup()

# Filtered to remove extended MHC (chr6, 25Mb to 34Mb).
  gene_coordinates <- 
  read_tsv("/share/pub/dengcy/Singlecell/COVID19/MAGMA/NCBI/NCBI37.3.gene.loc.extendedMHCexcluded",
           col_names = FALSE,col_types = 'cciicc') %>%
  mutate(start=ifelse(X3-100000<0,0,X3-100000),end=X4+100000) %>%
   select(X2,start,end,6,1) %>%
  rename(chr="X2", Gene="X6",ENTREZ="X1") %>%
  mutate(chr=paste0("chr",chr))

  exp2 <- inner_join(exp,gene_coordinates,by="Gene")

#Get number of genes that represent 10% of the dataset
	n_genes <- length(unique(exp2$ENTREZ))
	n_genes_to_keep <- (n_genes * 0.1) %>% round()

exp2 %>% filter(Expr_sum_mean>1) %>% ldsc_bedfile("column")
print("sucess!")
}

write_group  = function(df,Cell_type) {
  df <- select(df,column,chr,start,end,ENTREZ)
  dir.create(paste0("LDSC/Bed"), showWarnings = FALSE,recursive = TRUE)
  write_tsv(df[-1],paste0("LDSC/Bed/",make.names(unique(df[1])),".bed"),col_names = F)
return(df)
}

ldsc_bedfile <- function(d,Cell_type){
  d_spe <- d %>% group_by_(Cell_type) %>% top_n(.,n_genes_to_keep,specificity) 
  d_spe %>% do(write_group(.,Cell_type))
}


top10_function(liver_expr)
