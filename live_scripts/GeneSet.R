library(tidyverse)
library(clusterProfiler)

# over representation - KEGG

search_kegg_organism('mouse', by='common_name')

shrink.d11 <- readRDS("RObjects/Shrunk_Results.d11.rds")

sigGenes <- shrink.d11 %>%
  drop_na(Entrez, padj) %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  pull(Entrez)

keggRes <- enrichKEGG(gene = sigGenes, organism = 'mmu')

as_tibble(keggRes)

library(pathview)
logFC <- shrink.d11$log2FoldChange
names(logFC) <- shrink.d11$Entrez
pathview(gene.data = logFC,
         pathway.id = "mmu04612",
         species = "mmu",
         limit = list(gene = 20, cpd = 1))

# Exercise 1

#“mmu04659” or “mmu04658”

pathview(gene.data = logFC,
         pathway.id = "mmu04659",
         species = "mmu",
         limit = list(gene = 5, cpd = 1))

# ORA - GO terms

library(org.Mm.eg.db)

sigGenes_GO <- shrink.d11 %>%
  drop_na(padj) %>%
  filter(padj < 0.01 & abs(log2FoldChange) > 2) %>%
  pull(GeneID)

universe <- shrink.d11$GeneID

ego <- enrichGO(gene = sigGenes_GO,
                universe = universe,
                OrgDb = org.Mm.eg.db,
                keyType = "ENSEMBL",
                ont = "BP",
                pvalueCutoff = 0.01,
                readable = TRUE)

barplot(ego, showCategory = 20)
dotplot(ego, font.size =14)
library(enrichplot)
ego_pt <- pairwise_termsim(ego)
emapplot(ego_pt, category_label = 0.8)

# GSEA

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("msigdb")

library(msigdb)
library(ExperimentHub)

eh = ExperimentHub()
query(eh, c('msigdb', 'mm', '2023'))

msigdb.mm <- getMsigdb(org = 'mm', id = 'EZID', version = '2023.1')

listCollections(msigdb.mm)

hallmarks <- subsetCollection(msigdb.mm, 'h')
msigdb_ids <- geneIds(hallmarks)
term2gene <- enframe(msigdb_ids, name = "gs_name", value = "entrez") %>%
  unnest(entrez)

rankedGenes <- shrink.d11 %>%
  drop_na(GeneID, padj, log2FoldChange) %>%
  mutate(rank = log2FoldChange) %>%
  filter(!is.na(Entrez)) %>%
  arrange(desc(rank)) %>%
  pull(rank, Entrez)


gseaRes <- GSEA(rankedGenes,
                TERM2GENE = term2gene,
                pvalueCutoff = 1,
                minGSSize = 15,
                maxGSSize = 500)

View(as_tibble(gseaRes) %>%
       arrange(desc(abs(NES))) %>%
       dplyr::select(-core_enrichment) %>%
       mutate(across(c('enrichmentScore', 'NES'), ~round(.x, digit=3))) %>%
       mutate(across(c('pvalue', 'p.adjust', 'qvalue'), scales::scientific)))
gseaplot(gseaRes,
         geneSetID = 'HALLMARK_INFLAMMATORY_RESPONSE',
         title = 'HALLMARK_INFLAMMATORY_RESPONSE')

# Exercise 2

# part 1

rankedGenes.e11 <- shrink.d11 %>%
  drop_na(Entrez, pvalue, log2FoldChange) %>%
  mutate(rank = -log10(pvalue) * sign(log2FoldChange)) %>%
  arrange(desc(rank)) %>%
  pull(rank, Entrez)

# part 2

gseaRes.e11 <- GSEA(rankedGenes.e11,
                    TERM2GENE = term2gene,
                    pvalueCutoff = 1,
                    minGSSize = 15)

View(as_tibble(gseaRes.e11) %>%
       arrange(desc(abs(NES))) %>%
       dplyr::select(-core_enrichment) %>%
       mutate(across(c('enrichmentScore', 'NES'), ~round(.x, digit=3))) %>%
       mutate(across(c('pvalue', 'p.adjust', 'qvalue'), scales::scientific)))

# part 3

shrink.d33 <- readRDS('RObjects/Shrunk_Results.d33.rds')

rankedGenes.e33 <- shrink.d33 %>%
  drop_na(Entrez, pvalue, log2FoldChange) %>%
  mutate(rank = -log10(pvalue) * sign(log2FoldChange)) %>%
  arrange(desc(rank)) %>%
  pull(rank, Entrez)

gseaRes.e33 <- GSEA(rankedGenes.e33,
                    TERM2GENE = term2gene,
                    pvalueCutoff = 1,
                    minGSSize = 15)

tab.d33 <- as_tibble(gseaRes.e33) %>%
       arrange(desc(abs(NES))) %>%
       dplyr::select(-core_enrichment) %>%
       mutate(across(c('enrichmentScore', 'NES'), ~round(.x, digit=3))) %>%
       mutate(across(c('pvalue', 'p.adjust', 'qvalue'), scales::scientific))
