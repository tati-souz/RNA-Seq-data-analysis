
setwd("/Users/tatiane/Library/CloudStorage/OneDrive-UniversidadeFederaldoPará-UFPA/DOUTORADO/2.RESULTADOS_TRANSCRIPTOMA/Parkinson_Tati/Coexpressao_DifCoexp_dp/Results_DiffCoexp/")

########################################

# _ DCL _: Ligações Coexpressas Diferenciais - pares de genes com significância estatística em relação 
# à diferença dos coeficientes de correlação sob duas condições

# _ DCG _: Genes Coexpressos Diferenciais - genes com particularmente mais DCLs do que o esperado 
# por acaso

########################################

# CT vs CASE_DYSKINESIA
DCG_DPd <- read.delim("DCG_res_CT_DPd.txt", header = T, row.names = NULL)
DCL_DPd <- read.delim("DCL_res_CT_DPd.txt", header = T, row.names = NULL)

# CT vs CASE_NO_DYSKINESIA
DCG_DPnd <- read.delim("DCG_res_CT_DPnd.txt", header = T, row.names = NULL)
DCL_DPnd <- read.delim("DCL_res_CT_DPnd.txt", header = T, row.names = NULL)

########################################

# FUNÇÃO PARA ACHAR O SÍMBOLO DOS GENES

library(org.Hs.eg.db)
library(AnnotationDbi)

# Ensembl IDs To Gene Symbol
map_ensembl_to_symbol <- function(ensembl_ids) {
  gene_symbols <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = ensembl_ids, keytype = "ENSEMBL", column = "SYMBOL")
  return(gene_symbols)
}

##########################################
# Adicionar símbolos no DCLs DP-DYSKINESIA

# Search gene symbol to DCL_DPd
DCL_DPd

gene_symbol1 <- map_ensembl_to_symbol(DCL_DPd$Gene.1)
DCL_DPd$SymbolGene_1 <- gene_symbol1
gene_symbol2 <- map_ensembl_to_symbol(DCL_DPd$Gene.2)
DCL_DPd$SymbolGene_2 <- gene_symbol2
# Insert the ensembl value in gene_symbol1 when no matchs in map_ID
for (i in 1:nrow(DCL_DPd)) {
  # verificar se a linha está em branco
  if (is.na(DCL_DPd[i, 13]) || DCL_DPd[i,13] == "") {
    # adicionar o valor da coluna ao lado
    DCL_DPd[i,13] <- DCL_DPd[i,1]
  }
}
# Insert the ensembl value in gene_symbol1 when no matchs in map_ID
for (i in 1:nrow(DCL_DPd)) {
  # verificar se a linha está em branco
  if (is.na(DCL_DPd[i, 14]) || DCL_DPd[i,14] == "") {
    # adicionar o valor da coluna ao lado
    DCL_DPd[i,14] <- DCL_DPd[i,2]
  }
}

GeneList_DPd <- c(DCL_DPd$Gene.1, DCL_DPd$Gene.2)
write.table(GeneList_DPd, "./GeneList_DPd", quote = F, row.names = F)
write.csv(DCL_DPd, "/Users/tatiane/Library/CloudStorage/OneDrive-UniversidadeFederaldoPará-UFPA/DOUTORADO/2.RESULTADOS_TRANSCRIPTOMA/Parkinson_Tati/Coexpressao_DifCoexp_dp/Results_DiffCoexp/DCL_res_DPd_GeneSymbol.csv",
          quote = F, row.names = F)


###############################
# Adicionar símbolos no DCLs DP-DYSKINESIA

# Search gene symbol to DCL_DPd
DCL_DPnd

gene_symbol1 <- map_ensembl_to_symbol(DCL_DPnd$Gene.1)
DCL_DPnd$SymbolGene_1 <- gene_symbol1
gene_symbol2 <- map_ensembl_to_symbol(DCL_DPnd$Gene.2)
DCL_DPnd$SymbolGene_2 <- gene_symbol2
# Insert the ensembl value in gene_symbol1 when no matchs in map_ID
for (i in 1:nrow(DCL_DPnd)) {
  # verificar se a linha está em branco
  if (is.na(DCL_DPnd[i, 13]) || DCL_DPnd[i,13] == "") {
    # adicionar o valor da coluna ao lado
    DCL_DPnd[i,13] <- DCL_DPnd[i,1]
  }
}
# Insert the ensembl value in gene_symbol1 when no matchs in map_ID
for (i in 1:nrow(DCL_DPnd)) {
  # verificar se a linha está em branco
  if (is.na(DCL_DPnd[i, 14]) || DCL_DPnd[i,14] == "") {
    # adicionar o valor da coluna ao lado
    DCL_DPnd[i,14] <- DCL_DPnd[i,2]
  }
}

GeneList_DPnd <- c(DCL_DPnd$Gene.1, DCL_DPnd$Gene.2)
write.table(GeneList_DPnd, "./GeneList_DPnd", quote = F, row.names = F)
write.csv(DCL_DPnd, "/Users/tatiane/Library/CloudStorage/OneDrive-UniversidadeFederaldoPará-UFPA/DOUTORADO/2.RESULTADOS_TRANSCRIPTOMA/Parkinson_Tati/Coexpressao_DifCoexp_dp/Results_DiffCoexp/DCL_res_DPnd_GeneSymbol.csv",
          quote = F, row.names = F)





##############################################
# CHECANDO INTERAÇÕES ÚNICAS DENTRO DE CADA DCL

head(DCL_DPd)
head(DCL_DPnd)
library(dplyr)

### JEITO 1
# DPd e DPnd - Checando as interações EXCLUSIVAS (Independente da ordem) - Todas foram exclusivas
# Encontrar interações exclusivas em DCL_DPd
DCL_exclusivo_DPd <- anti_join(DCL_DPd, DCL_DPnd, by = c("Gene.1", "Gene.2")) %>%
  anti_join(DCL_DPnd, DCL_DPd, by = c("Gene.2", "Gene.1"))
# Encontrar interações exclusivas em DCL_DPnd
DCL_exclusivo_DPnd <- anti_join(DCL_DPnd, DCL_DPd, by = c("Gene.1", "Gene.2")) %>%
  anti_join(DCL_DPd, DCL_DPnd, by = c("Gene.2", "Gene.1"))

### JEITO 2
# Fazendo um upset das interações
library(UpSetR)

DCL_DPd <- DCL_DPd %>% 
  mutate(junto = apply(DCL_DPd[, 13:14], 1, paste, collapse = "/"))
DCL_DPd_list <- DCL_DPd$junto
DCL_DPnd <- DCL_DPnd %>% 
  mutate(junto = apply(DCL_DPnd[, 13:14], 1, paste, collapse = "/"))
DCL_DPnd_list <- DCL_DPnd$junto

listInput <- list("Interactions PD-D vs CT"= DCL_DPd_list,
                  "Interactions PD-ND vs CT"= DCL_DPnd_list)
up <- upset(fromList(listInput), order.by = "freq", matrix.color = "dodgerblue1", 
            nsets = 2, 
            point.size = 4,
            line.size = 1.5,
            sets.x.label = "Number of interactions",
            text.scale = c(1.5, 3, 1.5, 1.5, 3, 3))
up
library("ggvenn")
venn <- ggvenn(listInput,
                fill_color = c("deepskyblue", "dodgerblue1"),
                stroke_size = 0.8,
                set_name_size = 8,
                text_size = 16,
                show_percentage = F)
venn <- venn + ggtitle("Intersection of the DCLs") +
  theme(plot.title = element_text(size = 25, hjust = 0.5, face = "bold"))
venn




###################################################
### ANÁLISE DE SIMILARIDADE - Índice de Jaccard ###

# Função para calcular o índice de Jaccard
jaccard_index <- function(conjunto_A, conjunto_B) {
  # Transforma os conjuntos em vetores únicos
  vetor_A <- unique(conjunto_A)
  vetor_B <- unique(conjunto_B)
  # Tamanho da interseção
  intersecao <- length(intersect(vetor_A, vetor_B))
  # Tamanho da união
  uniao <- length(union(vetor_A, vetor_B))
  # Índice de Jaccard
  jaccard <- intersecao / uniao
  return(jaccard)
}

# Calcular o índice de Jaccard (INTERAÇÕES)
# DCL_DPd e DCL_DPnd
jaccard_interacoes <- jaccard_index(DCL_DPd$junto, DCL_DPnd$junto)
print(paste("O índice de Jaccard é:", jaccard_interacoes))
intersecao_interacoes <- intersect(DCL_DPd$junto, DCL_DPnd$junto)
intersecao_interacoes



###################################################
# FILTRANDO OS DCLS PELOS GENES DO DCG

### DP-D 
# Após filtrar por DCG, o DCL tem 23 interações (OU)
DCL_DPd_Filtrado <- DCL_DPd[DCL_DPd$Gene.1 %in% DCG_DPd$Gene |
                                 DCL_DPd$Gene.2 %in% DCG_DPd$Gene, ]
# 29 genes unicos
Genes_DCL_DPd_Filtrado <- c(DCL_DPd_Filtrado$Gene.1, DCL_DPd_Filtrado$Gene.2) %>% 
  unique() 
write.csv(DCL_DPd_Filtrado, file = "./DCL_DPd_Filtrado_por_DCG.csv", quote = F, row.names = F)

### DP-ND
# Após filtrar por DCG, o DCL tem 48 interações (E)
DCL_DPnd_Filtrado <- DCL_DPnd[DCL_DPnd$Gene.1 %in% DCG_DPnd$Gene |
                                 DCL_DPnd$Gene.2 %in% DCG_DPnd$Gene, ]
# 73 genes únicos
Genes_DCL_DPnd_Filtrado <- c(DCL_DPnd_Filtrado$Gene.1, DCL_DPnd_Filtrado$Gene.2) %>% 
  unique()
write.csv(DCL_DPnd_Filtrado, file = "./DCL_DPnd_Filtrado_por_DCG.csv", quote = F, row.names = F)


###################################################
# ANÁLISE DE ENRIQUECIMENTO DOS DCLs

### Verificando a lista de genes em cada DCL
length(GeneList_DPd) # 158 genes no DCL DP-D
length(GeneList_DPnd) # 6110 genes no DCL DP-ND

### Gene Ontology - Cluster Profiler
library(clusterProfiler)

############# DP-D #############
## SEA GO
go_enrich_DPd <- enrichGO(
  gene = GeneList_DPd,
  OrgDb = "org.Hs.eg.db",
  ont = "BP",
  pvalueCutoff = 0.5,
  qvalueCutoff = 0.5,
  pAdjustMethod = "BH",
  keyType = 'ENSEMBL',
  readable = T
)

# Obter dataframe
df_go_DPd <- go_enrich_DPd@result
df_go_DPd <- df_go_DPd[df_go_DPd$pvalue < 0.05,]
write.csv(df_go_DPd, "./GO_DCL_DPd.csv", quote = F, row.names = F)

############# DP-ND #############
## SEA GO
go_enrich_DPnd <- enrichGO(
  gene = GeneList_DPnd,
  OrgDb = "org.Hs.eg.db",
  ont = "BP",
  pvalueCutoff = 0.5,
  qvalueCutoff = 0.5,
  pAdjustMethod = "BH",
  keyType = 'ENSEMBL',
  readable = T
)

# Obter dataframe
df_go_DPnd <- go_enrich_DPnd@result
df_go_DPnd <- df_go_DPnd[df_go_DPnd$p.adjust < 0.05,]
write.csv(df_go_DPnd, "./GO_DCL_DPnd.csv", quote = F, row.names = F)




########### Intersection of pathways #############
nrow(df_go_DPd)
nrow(df_go_DPnd)
comuns <- intersect(df_go_DPd$ID, df_go_DPnd$ID)
length(comuns)


listInput <- list("Pathways PD-D vs CT"= df_go_DPd$ID,
                  "Pathways PD-ND vs CT"= df_go_DPnd$ID)
up <- upset(fromList(listInput), order.by = "freq", matrix.color = "dodgerblue1", 
            nsets = 2, 
            point.size = 4,
            line.size = 1.5,
            sets.x.label = "Number of interactions",
            text.scale = c(1.5, 3, 1.5, 1.5, 3, 3))
up

library(dplyr)
# Vias exclusivas
# DP-D - 110 vias
ids_exclusivos_DPd <- setdiff(df_go_DPd$ID, comuns)
df_go_DPd_exclusivos <- df_go_DPd %>% filter(ID %in% ids_exclusivos_DPd)
write.csv(df_go_DPd_exclusivos, "./GO_DCL_Exclusivos_DPd.csv", quote = F, row.names = F)

# DP-ND - 178 vias
ids_exclusivos_DPnd <- setdiff(df_go_DPnd$ID, comuns)
df_go_DPnd_exclusivos <- df_go_DPnd %>% filter(ID %in% ids_exclusivos_DPnd)
write.csv(df_go_DPnd_exclusivos, "./GO_DCL_Exclusivos_DPnd.csv", quote = F, row.names = F)



###############################
#### plot vias exclusivas #####

##############
library(multienrichjam)

# DP-D
go_enrich_DPd <- enrichDF2enrichResult(df_go_DPd_exclusivos)

# Barplot - GO
bar <- barplot(go_enrich_DPd,
               drop = TRUE, 
               showCategory = 10, 
               title = "Biological Processes in the PD-D group",
               font.size = 14)
bar
ggsave("./GO_DPd_Exclusive_DCL_Barplot.pdf", plot = bar, 
       width = 8, height = 10, dpi = 600)

# Dotplot - GO
dot <- dotplot(go_enrich_DPd, showCategory = 10,  
               title = "Biological Processes - DPd group")
dot
ggsave("./GO_DPd_Exclusive_DCL_Dotplot.png", plot = dot, 
       width = 8, height = 6, dpi = 600)

# Cnetplot - GO
cnet <- cnetplot(go_enrich_DPd, title = "Biological Processes - DPd group",
                 cex.params = list(category_label = 0.9, gene_label = 0.8),
                 color.params = list(gene = "cornflowerblue", category = "green3"),
                 showCategory = 10,
                 layout = "kk",
                 circular = F,
                 colorEdge = F)
cnet
ggsave("./GO_DPd_Exclusive_DCL_Cnetplot.png", plot = cnet, 
       width = 8, height = 6, dpi = 600)

# Emapplot
library(enrichplot)
x <- pairwise_termsim(go_enrich_DPd)
emap <- emapplot(x, layout = "kk",cex.params = list(category_label = 0.7), 
                 title = "GO Biological Pathways - Emapplot")
emap
ggsave("GO_DPd_Exclusive_DCL_Emapplot.png", plot = emap, 
       width = 8, height = 6, dpi = 600)
# Heatplot
heatplot(go_enrich_DPd, showCategory = 10)
ggsave("GO_DPd_Exclusive_DCL_heatplot.png", plot = cnet, 
       width = 8, height = 6, dpi = 600)


##############

# DP-ND
go_enrich_DPnd <- enrichDF2enrichResult(df_go_DPnd_exclusivos)

# Barplot - GO
bar <- barplot(go_enrich_DPnd,
               drop = TRUE, 
               showCategory = 10, 
               title = "Biological Processes in the PD-ND group")
bar
ggsave("./GO_DPnd_Exclusive_DCL_Barplot.png", plot = bar, 
       width = 8, height = 6, dpi = 600)

# Dotplot - GO
dot <- dotplot(go_enrich_DPnd, showCategory = 10,  
               title = "Biological Processes - DPnd group")
dot
ggsave("./GO_DPnd_Exclusive_DCL_Dotplot.png", plot = dot, 
       width = 8, height = 6, dpi = 600)

# Cnetplot - GO
cnet <- cnetplot(go_enrich_DPnd, title = "Biological Processes - DPnd group",
                 cex.params = list(category_label = 0.6, gene_label = 0.6),
                 color.params = list(gene = "cornflowerblue", category = "green3"),
                 showCategory = 10,
                 layout = "kk",
                 circular = F,
                 colorEdge = F)
cnet
ggsave("./GO_DPnd_Exclusive_DCL_Cnetplot.png", plot = cnet, 
       width = 8, height = 6, dpi = 600)


