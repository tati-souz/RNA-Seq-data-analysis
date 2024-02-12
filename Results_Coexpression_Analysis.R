########################################################
# ANÁLISE DOS RESULTADOS DA COEXPRESSÃO
########################################################

setwd("/Users/tatiane/Library/CloudStorage/OneDrive-UniversidadeFederaldoPará-UFPA/DOUTORADO/2.RESULTADOS_TRANSCRIPTOMA/Parkinson_Tati/Coexpressao_DifCoexp_dp/Resultados_Coexpressão")

##########################################
# ANÁLISE DE ENRIQUECIMENTO ##############
##########################################
library(org.Hs.eg.db)
### GRUPO DP-COM-DISCINESIA

########## ORA OF INTERACTION GENES 
ora_discinesia <- read.delim("./Cemitool_Caso_Discinesia/Tables/ora.tsv")
ora_discinesia <- ora_discinesia[ora_discinesia$p.adjust < 0.05,]

ora_discinesia_signif <- ora_discinesia[ora_discinesia$Module %in% c("M4","M5","M6"),]

      # Função para mapear Ensembl para símbolos de gene
      map_ensembl_to_symbol <- function(ensembl) {
        gene_symbols <- mapIds(org.Hs.eg.db, keys = ensembl, keytype = "ENSEMBL", column = "SYMBOL")
        return(gene_symbols)
      }
      # Aplicar a função de mapeamento para cada linha da tabela de resultados
      ora_discinesia_signif$SymbolGenes <- sapply(strsplit(ora_discinesia_signif$geneID, "/"), map_ensembl_to_symbol)
      ora_discinesia_signif$SymbolGenes <- gsub(", ", "/", ora_discinesia_signif$SymbolGenes)
      
      # Usando regex para extrair apenas os símbolos de gene
      symbol_genes <- sapply(strsplit(ora_discinesia_signif$SymbolGenes, "/"), function(pair) {
        symbols <- sapply(pair, function(item) {
          symbol <- gsub("[0-9]+ = \"([A-Z0-9]+)\"", "\\1", item)
          return(symbol)
        })
        symbols <- symbols[symbols != ""]
        return(paste(symbols, collapse = "/"))
      })
      # Substituir a coluna SymbolGenes pelos símbolos extraídos
      ora_discinesia_signif$SymbolGenes <- symbol_genes
      ora_discinesia_signif$SymbolGenes <- sub("c\\(|\\)", "", ora_discinesia_signif$SymbolGenes)
      ora_discinesia_signif$SymbolGenes <- gsub("(?<=/|^)ENSG", "", ora_discinesia_signif$SymbolGenes, perl = TRUE)
      ora_discinesia_signif$SymbolGenes <- gsub("\\)", "", ora_discinesia_signif$SymbolGenes, perl = TRUE)
      
write.table(ora_discinesia_signif, "./Cemitool_Caso_Discinesia/Tables/ora_genesymbol.txt")


### GRUPO DP-SEM-DISCINESIA

########## ORA OF INTERACTION GENES 
ora_no_discinesia <- read.delim("./Cemitool_Caso_Sem_Discinesia/Tables/ora.tsv")
ora_no_discinesia <- ora_no_discinesia[ora_no_discinesia$p.adjust < 0.05,]

ora_no_discinesia_signif <- ora_no_discinesia[ora_no_discinesia$Module == "M2",]

# Aplicar a função de mapeamento para cada linha da tabela de resultados (Função criada acima)
ora_no_discinesia_signif$SymbolGenes <- sapply(strsplit(ora_no_discinesia_signif$geneID, "/"), map_ensembl_to_symbol)
ora_no_discinesia_signif$SymbolGenes <- gsub(", ", "/", ora_no_discinesia_signif$SymbolGenes)

# Usando regex para extrair apenas os símbolos de gene
symbol_genes <- sapply(strsplit(ora_no_discinesia_signif$SymbolGenes, "/"), function(pair) {
  symbols <- sapply(pair, function(item) {
    symbol <- gsub("[0-9]+ = \"([A-Z0-9]+)\"", "\\1", item)
    return(symbol)
  })
  symbols <- symbols[symbols != ""]
  return(paste(symbols, collapse = "/"))
})
# Substituir a coluna SymbolGenes pelos símbolos extraídos
ora_no_discinesia_signif$SymbolGenes <- symbol_genes
ora_no_discinesia_signif$SymbolGenes <- sub("c\\(|\\)", "", ora_no_discinesia_signif$SymbolGenes)
ora_no_discinesia_signif$SymbolGenes <- gsub("(?<=/|^)ENSG", "", ora_no_discinesia_signif$SymbolGenes, perl = TRUE)
ora_no_discinesia_signif$SymbolGenes <- gsub("\\)", "", ora_no_discinesia_signif$SymbolGenes, perl = TRUE)

write.table(ora_no_discinesia_signif, "./Cemitool_Caso_Sem_Discinesia/Tables/ora_genesymbol.txt")



### GRUPO CT

########## ORA OF INTERACTION GENES 
ora_controle <- read.delim("./Cemitool_Controle/Tables/ora.tsv")
ora_controle <- ora_controle[ora_controle$p.adjust < 0.05,]

ora_controle_signif <- ora_controle[ora_controle$Module == "M2",]

# Aplicar a função de mapeamento para cada linha da tabela de resultados (Função criada acima)
ora_controle_signif$SymbolGenes <- sapply(strsplit(ora_controle_signif$geneID, "/"), map_ensembl_to_symbol)
ora_controle_signif$SymbolGenes <- gsub(", ", "/", ora_controle_signif$SymbolGenes)

# Usando regex para extrair apenas os símbolos de gene
symbol_genes <- sapply(strsplit(ora_controle_signif$SymbolGenes, "/"), function(pair) {
  symbols <- sapply(pair, function(item) {
    symbol <- gsub("[0-9]+ = \"([A-Z0-9]+)\"", "\\1", item)
    return(symbol)
  })
  symbols <- symbols[symbols != ""]
  return(paste(symbols, collapse = "/"))
})
# Substituir a coluna SymbolGenes pelos símbolos extraídos
ora_controle_signif$SymbolGenes <- symbol_genes
ora_controle_signif$SymbolGenes <- sub("c\\(|\\)", "", ora_controle_signif$SymbolGenes)
ora_controle_signif$SymbolGenes <- gsub("(?<=/|^)ENSG", "", ora_controle_signif$SymbolGenes, perl = TRUE)
ora_controle_signif$SymbolGenes <- gsub("\\)", "", ora_controle_signif$SymbolGenes, perl = TRUE)

write.table(ora_controle_signif, "./Cemitool_Controle/Tables/ora_genesymbol.txt")



### GRUPO GERAL

########## ORA OF INTERACTION GENES
ora_geral <- read.delim("./Cemitool_Geral/Tables/ora.tsv")
ora_geral <- ora_geral[ora_geral$p.adjust < 0.05,]

ora_geral_signif <- ora_geral[ora_geral$Module %in% c("M1","M2"),]

# Aplicar a função de mapeamento para cada linha da tabela de resultados
ora_geral_signif$SymbolGenes <- sapply(strsplit(ora_geral_signif$geneID, "/"), map_ensembl_to_symbol)
ora_geral_signif$SymbolGenes <- gsub(", ", "/", ora_geral_signif$SymbolGenes)

# Usando regex para extrair apenas os símbolos de gene
symbol_genes <- sapply(strsplit(ora_geral_signif$SymbolGenes, "/"), function(pair) {
  symbols <- sapply(pair, function(item) {
    symbol <- gsub("[0-9]+ = \"([A-Z0-9]+)\"", "\\1", item)
    return(symbol)
  })
  symbols <- symbols[symbols != ""]
  return(paste(symbols, collapse = "/"))
})
# Substituir a coluna SymbolGenes pelos símbolos extraídos
ora_geral_signif$SymbolGenes <- symbol_genes
ora_geral_signif$SymbolGenes <- sub("c\\(|\\)", "", ora_geral_signif$SymbolGenes)
ora_geral_signif$SymbolGenes <- gsub("(?<=/|^)ENSG", "", ora_geral_signif$SymbolGenes, perl = TRUE)
ora_geral_signif$SymbolGenes <- gsub("\\)", "", ora_geral_signif$SymbolGenes, perl = TRUE)

write.table(ora_geral_signif, "./Cemitool_Geral/Tables/ora_genesymbol.txt")


#######################################
# ANÁLISE VIAS EXCLUSIVAS: DP-D E DP-ND
#######################################
# Não houve via exclusiva de DP-ND, apenas para DP-D (24 vias).
# Das 24 vias exclusivas de DP-D, 10 possielmente relacionadas com o fenótipo foram selecionadas discussão

### verificando as interseções das vias
library(UpSetR)
listInput1 <- list("Pathways PD-D"= ora_discinesia_signif$Description,
                   "Pathways PD-ND"= ora_no_discinesia_signif$Description,
                   "Pathways CT"= ora_controle_signif$Description,
                   "Pathways CT + PD-D + PD-ND"= ora_geral_signif$Description)

up1 <- upset(fromList(listInput1[c(1,2,3,4)]), order.by = "freq", matrix.color = "dodgerblue1", 
             nsets = 4, 
             point.size = 4,
             line.size = 1.5,
             sets.x.label = "Number of pathways",
             text.scale = c(1.5, 3, 1.5, 1.5, 3, 3))
up1

# 24 vias exclusivas na DP-D
exclusive_path_DPd <- setdiff(ora_discinesia_signif$Description, ora_geral_signif$Description)
exclusive_path_DPd <- setdiff(exclusive_path_DPd, ora_controle_signif$Description)
exclusive_path_DPd <- setdiff(exclusive_path_DPd, ora_no_discinesia_signif$Description)
length(exclusive_path_DPd)

ora_exclusive_discinesia_signif <- ora_discinesia_signif[ora_discinesia_signif$Description %in% exclusive_path_DPd,-11]
write.csv(ora_exclusive_discinesia_signif, "./ora_ExclusivePathways_DP-D.csv")

# Barplot das 24 vias exclusivas - PD-D
# Transformar ora_exclusive_discinesia_signif em objeto de enriquecimento 
library(multienrichjam)
ora_exc_DPd <- ora_exclusive_discinesia_signif
colnames(ora_exc_DPd)[c(9,11)] <- c("ensengID", "geneID")
ora_exc_DPd_obj <- enrichDF2enrichResult(ora_exc_DPd)

bar <- barplot(ora_exc_DPd_obj,
               drop = TRUE, 
               showCategory = 24, 
               title = "ORA - Exclusive DP-D pathways",
               font.size = 9)
bar
library(ggplot2)
ggsave("./Ora_ExclusivePathways_DP-D_Barplot.png", plot = bar, 
       width = 8, height = 12, dpi = 600)

# Selecionando 10 vias das 24 que eu pretendo discutir:
vias_selecionadas <- c(" Metabolism of amino acids and derivatives",
                       " Axon guidance",
                       " Nervous system development",
                       " Cellular responses to stress",
                       " L1CAM interactions",
                       " Extracellular matrix organization")

ora_exclusive_discinesia_discussion <- ora_exclusive_discinesia_signif[ora_exclusive_discinesia_signif$ID %in% vias_selecionadas,]
write.csv(ora_exclusive_discinesia_discussion, "./ora_ExclusivePathways_DP-D_Discussion.csv")

# Transformar ora_exclusive_discinesia_discussion em objeto de enriquecimento 
library(multienrichjam)
ora_exc_DPd_select <- ora_exclusive_discinesia_discussion
ora_exc_DPd_select_obj <- enrichDF2enrichResult(ora_exc_DPd_select)

# PLOTS

# Barplot - ORA
bar <- barplot(ora_exc_DPd_select_obj,
               drop = TRUE, 
               showCategory = 10, 
               title = "ORA - key exclusive pathways in PD-D")
bar
ggsave("./Ora_ExclusivePathways_DP-D_Discussion_Barplot.png", plot = bar, 
       width = 8, height = 6, dpi = 600)

# Cnetplot
cnet <- cnetplot(ora_exc_DPd_select_obj, title = "ORA Dyskinesia - Cnetplot", showCategory = 10,
                 cex.params = list(category_label = 0.6, gene_label = 0.6),
                 layout = "fr",
                 circular = F,
                 color.params = list(edge = F, gene = "cornflowerblue", category = "green3"))
cnet
ggsave("./Ora_ExclusivePathways_DP-D_Discussion_Cnetplot.png", plot = cnet, 
       width = 8, height = 6, dpi = 600)





###############################################################################
###############################################################################
###############################################################################

# ANÁLISE DAS REDES DE CO-EXPRESSÃO

###################
# INTERACTIONS DP-D
interactions_discinesia <- read.delim("./Cemitool_Caso_Discinesia/Tables/interactions.tsv")

# Função para mapear Entrez IDs para símbolos de gene
library(clusterProfiler)
library(org.Hs.eg.db)
map_ensembl_to_symbol <- function(ensembl_ids) {
  gene_symbols <- mapIds(org.Hs.eg.db, keys = ensembl_ids, keytype = "ENSEMBL", column = "SYMBOL")
  return(gene_symbols)
}
# Aplicar a função de mapeamento para cada linha da tabela de resultados
interactions_discinesia_1 <- interactions_discinesia
interactions_discinesia_1$SymbolGene_1 <- sapply(interactions_discinesia_1$Gene1, map_ensembl_to_symbol)
interactions_discinesia_1$SymbolGene_2 <- sapply(interactions_discinesia_1$Gene2, map_ensembl_to_symbol)

write.csv(interactions_discinesia_1, "./Cemitool_Caso_Discinesia/Tables/interactions_GeneSymbol.tsv",quote = F, row.names = F)

###################
# INTERACTIONS DP-ND
interactions_sem_discinesia <- read.delim("./Cemitool_Caso_Sem_Discinesia/Tables/interactions.tsv")

# Aplicar a função de mapeamento para cada linha da tabela de resultados (Função definida acima)
interactions_sem_discinesia_1 <- interactions_sem_discinesia
interactions_sem_discinesia_1$SymbolGene_1 <- sapply(interactions_sem_discinesia_1$Gene1, map_ensembl_to_symbol)
interactions_sem_discinesia_1$SymbolGene_2 <- sapply(interactions_sem_discinesia_1$Gene2, map_ensembl_to_symbol)

write.csv(interactions_sem_discinesia_1, "./Cemitool_Caso_Sem_Discinesia/Tables/interactions_GeneSymbol.tsv",quote = F, row.names = F)

###################
# INTERACTIONS CT
interactions_controle <- read.delim("./Cemitool_Controle/Tables/interactions.tsv")

# Aplicar a função de mapeamento para cada linha da tabela de resultados (Função criada acima)
interactions_controle_1 <- interactions_controle
interactions_controle_1$SymbolGene_1 <- sapply(interactions_controle_1$Gene1, map_ensembl_to_symbol)
interactions_controle_1$SymbolGene_2 <- sapply(interactions_controle_1$Gene2, map_ensembl_to_symbol)

write.csv(interactions_controle_1, "./Cemitool_Controle/Tables/interactions_GeneSymbol.tsv",quote = F, row.names = F)

###################
# INTERACTIONS GERAL
interactions_geral <- read.delim("./Cemitool_Geral/Tables/interactions.tsv")

# Aplicar a função de mapeamento para cada linha da tabela de resultados
interactions_geral_1 <- interactions_geral
interactions_geral_1$SymbolGene_1 <- sapply(interactions_geral_1$Gene1, map_ensembl_to_symbol)
interactions_geral_1$SymbolGene_2 <- sapply(interactions_geral_1$Gene2, map_ensembl_to_symbol)

write.csv(interactions_geral_1, "./Cemitool_Geral/Tables/interactions_GeneSymbol.tsv", quote = F, row.names = F)


#### VERIFICANDO A INTERSEÇÃO DAS INTERAÇÕES ENTRE OS GRUPOS ANALISADOS
# Testei com a ordem inversa das interações e o resultado não mudou

library(dplyr)

DPd_interactions <- interactions_discinesia_1 %>% 
  mutate(junto = apply(interactions_discinesia_1[, 4:5], 1, paste, collapse = "/"))
DPd_interactions_list <- DPd_interactions$junto

DPnd_interactions <- interactions_sem_discinesia_1 %>% 
  mutate(junto = apply(interactions_sem_discinesia_1[, 4:5], 1, paste, collapse = "/"))
DPnd_interactions_list <- DPnd_interactions$junto

CT_interactions <- interactions_controle_1 %>% 
  mutate(junto = apply(interactions_controle_1[, 4:5], 1, paste, collapse = "/"))
CT_interactions_list <- CT_interactions$junto

Geral_interactions <- interactions_geral_1 %>% 
  mutate(junto = apply(interactions_geral_1[, 4:5], 1, paste, collapse = "/"))
Geral_interactions_list <- Geral_interactions$junto


library(UpSetR)
listInput <- list("Interactions PD-D"= DPd_interactions_list,
                  "Interactions PD-ND"= DPnd_interactions_list,
                  "Interactions CT"= CT_interactions_list,
                  "Interactions CT + PD-D + PD-ND"= Geral_interactions_list)

up <- upset(fromList(listInput), order.by = "freq", matrix.color = "dodgerblue1", 
            nsets = 4, 
            point.size = 4,
            line.size = 1.5,
            sets.x.label = "Number of parwaise interactions",
            text.scale = c(1.5, 3, 1.5, 1.5, 3, 3))
up

# Interações exclusivas na DP-D (12)
exclusive_inter_DPd <- setdiff(DPd_interactions$junto, DPnd_interactions$junto)
exclusive_inter_DPd <- setdiff(exclusive_inter_DPd, CT_interactions$junto)
exclusive_inter_DPd <- setdiff(exclusive_inter_DPd, Geral_interactions$junto)
length(exclusive_inter_DPd)

DPd_interactions_exclusive <- DPd_interactions[DPd_interactions$junto %in% exclusive_inter_DPd,]
DPd_interactions_exclusive <- DPd_interactions_exclusive[,4:5]
DPd_interactions_exclusive <- unique(DPd_interactions_exclusive)
Genelist_DPd_Coexp <- c(DPd_interactions_exclusive$SymbolGene_1, DPd_interactions_exclusive$SymbolGene_2) %>% unique()
length(Genelist_DPd_Coexp)

write.csv(DPd_interactions_exclusive, "./Cemitool_Caso_Discinesia/Tables/interactions_exclusive_GeneSymbol.csv", quote = F, row.names = F)
write.csv(Genelist_DPd_Coexp, "./Cemitool_Caso_Discinesia/Tables/Genelist_DPd_Coexp.csv", quote = F, row.names = F)


# Interações exclusivas na DP-ND (23)
exclusive_inter_DPnd <- setdiff(DPnd_interactions$junto, DPd_interactions$junto)
exclusive_inter_DPnd <- setdiff(exclusive_inter_DPnd, CT_interactions$junto)
exclusive_inter_DPnd <- setdiff(exclusive_inter_DPnd, Geral_interactions$junto)
length(exclusive_inter_DPnd)

DPnd_interactions_exclusive <- DPnd_interactions[DPnd_interactions$junto %in% exclusive_inter_DPnd,]
DPnd_interactions_exclusive <- DPnd_interactions_exclusive[,4:5]
DPnd_interactions_exclusive <- unique(DPnd_interactions_exclusive)
Genelist_DPnd_Coexp <- c(DPnd_interactions_exclusive$SymbolGene_1, DPnd_interactions_exclusive$SymbolGene_2) %>% unique()
length(Genelist_DPnd_Coexp)

write.csv(DPnd_interactions_exclusive, "./Cemitool_Caso_Sem_Discinesia/Tables/interactions_exclusive_GeneSymbol.csv", quote = F, row.names = F)
write.csv(Genelist_DPnd_Coexp, "./Cemitool_Caso_Sem_Discinesia/Tables/Genelist_DPnd_Coexp.csv", quote = F, row.names = F)





######### Carregando tabelas salvas
path1 <- read.csv("ora_ExclusivePathways_DP-D.csv", header = T)
path <- read.csv("ora_ExclusivePathways_DP-D_Discussion.csv", header = T)
# Função para mapear Ensembl para símbolos de gene
map_ensembl_to_symbol <- function(ensembl) {
  gene_symbols <- mapIds(org.Hs.eg.db, keys = ensembl, keytype = "ENSEMBL", column = "SYMBOL")
  return(gene_symbols)
}
# Aplicar a função de mapeamento para cada linha da tabela de resultados
path$SymbolGenes <- sapply(strsplit(path$geneID, "/"), map_ensembl_to_symbol)
path$SymbolGenes <- gsub(", ", "/", path$SymbolGenes)

# Usando regex para extrair apenas os símbolos de gene
symbol_genes <- sapply(strsplit(path$SymbolGenes, "/"), function(pair) {
  symbols <- sapply(pair, function(item) {
    symbol <- gsub("[0-9]+ = \"([A-Z0-9]+)\"", "\\1", item)
    return(symbol)
  })
  symbols <- symbols[symbols != ""]
  return(paste(symbols, collapse = "/"))
})
# Substituir a coluna SymbolGenes pelos símbolos extraídos
path$SymbolGenes <- symbol_genes
path$SymbolGenes <- sub("c\\(|\\)", "", path$SymbolGenes)
path$SymbolGenes <- gsub("(?<=/|^)ENSG", "", path$SymbolGenes, perl = TRUE)
path$SymbolGenes <- gsub("\\)", "", path$SymbolGenes, perl = TRUE)

metab_aminoacids <- path$SymbolGenes[1]
genes_metab_aminoacids <- strsp