if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq",force = TRUE)
BiocManager::install("dada2",force = TRUE)
library(readr);library(phyloseq);library("dada2");library(tidyverse);library(data.table);library(vegan);
library(Biostrings)

asv_tab <- read.table("~/Documents/asv_table.tsv", header=TRUE, row.names=1, sep="\t", check.names = FALSE)
asv_tab_t <- t(asv_tab)
asv_tab_t <- as.data.frame(asv_tab_t)
asv_tab_t$ASV <- rownames(asv_tab_t)

asv_seqs <- DNAStringSet(rownames(asv_tab_t))

taxa <- assignTaxonomy(asv_seqs, "./silva_nr99_v138.2_toGenus_trainset.fa", multithread=TRUE)
taxa_df <- as.data.frame(taxa)
taxa_df$ASV <- rownames(taxa_df)
taxa_df <- taxa_df %>% filter(Kingdom == "Bacteria")

df <- left_join(taxa_df, asv_tab_t, by="ASV")

df$Genus[is.na(df$Genus) | df$Genus == ""] <- "unclassified"

df_genus <- df %>%
  select(Genus, starts_with("ERR"))

asv.genus <- df_genus %>%
  group_by(Genus) %>%
  summarise(across(starts_with("ERR"), sum)) %>%
  ungroup()

asv.genus.sum <- asv.genus %>%
  group_by(Genus) %>%
  summarise(across(starts_with("ERR"), sum), .groups = "drop")

genus_mat <- asv.genus.sum %>%
  column_to_rownames("Genus") %>%
  t()

genus_ra <- genus_mat %>%
  apply(1, function(x) x / sum(x)) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id")

meta <- read.delim("~/Documents/accession_arm.tsv")
colnames(meta) <- c("sample_id","diagnosis")

plot_df <- genus_ra %>%
  left_join(meta, by = "sample_id")

longgenus <- plot_df %>%
  pivot_longer(cols = -c(sample_id, diagnosis),
               names_to = "Genus",
               values_to = "RA")

#####creating the phylogenetic object
otu_mat <- asv_tab_t %>%
  select(starts_with("ERR")) %>%        
  as.matrix()

### taxa matrix
tax_mat <- taxa_df %>%
  select(-ASV) %>%                       
  as.matrix()
rownames(tax_mat) <- taxa_df$ASV        

###new meta data for phylo seq
meta_ps <- meta %>%
  filter(sample_id %in% colnames(otu_mat))

rownames(meta_ps) <- meta_ps$sample_id

#####phyloseq object
ps <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(meta_ps)
)

#### alpha+richness and beta diversity
alpha_df <- estimate_richness(
  ps,
  measures = c("Shannon", "InvSimpson", "Chao1", "ACE")
)

alpha_df$sample_id <- rownames(alpha_df)

alpha_df <- alpha_df %>%
  left_join(meta_ps, by = "sample_id")

##### to check normality 
library(rstatix)
library(ggplot2)
library(ggpubr)
normality_results <- alpha_df %>%
  group_by(diagnosis) %>%
  shapiro_test(Shannon) #### not normal

p<- alpha_df %>%
  ggplot(aes(x = diagnosis, y = Shannon, fill = diagnosis)) +
  geom_boxplot() +
  theme_bw() +
  labs(title = "Shannon Diversity by Diagnosis")

stat_test_results <- alpha_df %>%
  wilcox_test(Shannon ~ diagnosis) %>%
  add_xy_position(x = "diagnosis")

# 
sigplotalpha<- p + stat_pvalue_manual(
  stat_test_results,
  label = "p", 
  tip.length = 0.01,  
  bracket.nudge.y = 0.5 
)
sigplotalpha


##### phylogenic diversity 
# BiocManager::install("DECIPHER")
library(DECIPHER)
library(phangorn)
library(picante)
library(ape)

seqs <- DNAStringSet(rownames(otu_mat))

alignment <- AlignSeqs(seqs)

phang.align <- phyDat(as.matrix(alignment), type = "DNA")

dm <- dist.ml(phang.align)

treeNJ <- NJ(dm)
treeNJ <- ladderize(treeNJ)

treeNJ$tip.label <- rownames(otu_mat)

ps_tree <- merge_phyloseq(ps, treeNJ)

otu_for_pd <- as(otu_table(ps_tree), "matrix")


if (taxa_are_rows(ps_tree)) {
  otu_for_pd <- t(otu_for_pd)
}

pd_res <- pd(otu_for_pd, phy_tree(ps_tree), include.root = FALSE)
pd_res$sample_id <- rownames(pd_res)

pd_df <- pd_res %>%
  left_join(meta_ps, by = "sample_id")


test_res <- wilcox.test(PD ~ diagnosis, data = pd_df)
p_val <- signif(test_res$p.value, 3)
p_val

ggplot(pd_df, aes(x = diagnosis, y = PD, fill = diagnosis)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
  annotate("text",
           x = 1.5,
           y = max(pd_df$PD) * 1.05,
           label = paste0("Wilcoxon p = ", p_val),
           size = 5) +
  labs(
    title = "Phylogenetic Diversity Across Groups",
    y = "Faith's PD",
    x = ""
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"))

## PERMANOVA and beta diversity
bray_dist <- vegan::vegdist(t(otu_table(ps)), method = "bray")
set.seed(123)
permanova_res <- adonis2(bray_dist ~ diagnosis, data = meta_ps, permutations = 999)
permanova_res
ord_pcoa <- ordinate(ps, method = "PCoA", distance = "bray")
###pcoa
bray_dist <- vegan::vegdist(t(otu_table(ps)), method = "bray")
pcoa_bc <- ape::pcoa(bray_dist)
pcoa_df <- as.data.frame(pcoa_bc$vectors[, 1:2])
pcoa_df$sample_id <- rownames(pcoa_df)

pcoa_df <- pcoa_df %>%
  left_join(meta_ps, by = "sample_id")
ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, color = diagnosis)) +
  geom_point(size = 4, alpha = 0.9) +
  theme_bw(base_size = 14) +
  labs(
    title = "PCoA showing Bray–Curtis Dissimilarity of AD Group and Control Group",
    x = paste0("PCoA 1 (", round(pcoa_bc$values$Relative_eig[1] * 100, 1), "%)"),
    y = paste0("PCoA 2 (", round(pcoa_bc$values$Relative_eig[2] * 100, 1), "%)")
  )
#### trying again but with centroids and ellipses
centroids <- pcoa_df %>%
  group_by(diagnosis) %>%
  summarise(
    Axis.1 = mean(Axis.1),
    Axis.2 = mean(Axis.2)
  )
ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, color = diagnosis)) +
  geom_point(size = 4, alpha = 0.9) +
  stat_ellipse(type = "t", linewidth = 1.2) +
  geom_point(
    data = centroids,
    aes(x = Axis.1, y = Axis.2, fill = diagnosis),
    color = "black",
    shape = 24,   
    size = 5
  ) +
  theme_bw(base_size = 14) +
  labs(
    title = "PCoA showing Bray–Curtis Dissimilarity of AD Group and Control Group",
    x = paste0("PCoA 1 (", round(pcoa_bc$values$Relative_eig[1] * 100, 1), "%)"),
    y = paste0("PCoA 2 (", round(pcoa_bc$values$Relative_eig[2] * 100, 1), "%)")
  ) 

### differential abundance analysis
genus_mat <- asv.genus.sum %>%
  column_to_rownames("Genus") %>%
  as.matrix()

genus_mat <- genus_mat[, meta_ps$sample_id]

ps_genus <- phyloseq(
  otu_table(genus_mat, taxa_are_rows = TRUE),
  sample_data(meta_ps)
)
#BiocManager::install("DESeq2")
library(DESeq2)
meta_ps$diagnosis <- factor(meta_ps$diagnosis)
levels(meta_ps$diagnosis)

ps_genus <- phyloseq(
  otu_table(genus_mat, taxa_are_rows = TRUE),
  sample_data(meta_ps)
)
dds <- phyloseq_to_deseq2(ps_genus, ~ diagnosis)
dds <- DESeq(dds, fitType = "local")
res <- results(dds, contrast = c("diagnosis", 
                                 "stool sample Control", 
                                 "stool sample AD"))
res_df <- as.data.frame(res)
res_df$Genus <- rownames(res_df)
sig_genera <- res_df %>%
  filter(padj < 0.05)
sig_genera

sig_list <- sig_genera$Genus
plot_df_sig <- longgenus %>%
  filter(Genus %in% sig_list)

##### corrected
ggplot(plot_df_sig, aes(x = diagnosis, y = RA, fill = diagnosis)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~ Genus, scales = "free_y") +
  theme_bw(base_size = 14) +
  ylab("Relative Abundance") +
  xlab("") +
  ggtitle("Differentially Abundant Genera (AD vs Control)") +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

sig_genera %>%
  mutate(direction = ifelse(log2FoldChange > 0, "Increased_in_AD", "Decreased_in_AD")) %>%
  arrange(padj)
####trying barplot
res_df_sorted <- res_df %>%
  arrange(log2FoldChange)
ggplot(res_df_sorted,
       aes(x = reorder(Genus, log2FoldChange),
           y = log2FoldChange,
           fill = log2FoldChange > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c("TRUE" = "firebrick",   
               "FALSE" = "steelblue")  
  ) +
  theme_bw(base_size = 14) +
  labs(
    title = "Log2 Fold Change for All Genera (AD vs Control)",
    x = "Genus",
    y = "Log2 Fold Change",
    fill = "Higher in AD?"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )
#####trying for uncorrected
genus_tests <- longgenus %>%
  group_by(Genus) %>%
  wilcox_test(RA ~ diagnosis) %>%
  select(Genus, p)

plot_df_uncorrected <- longgenus %>%
  left_join(genus_tests, by = "Genus") %>%   
  left_join(phylum_map,  by = "Genus")

direction_map <- res_df %>%
  select(Genus, log2FoldChange) %>%
  mutate(direction = ifelse(log2FoldChange > 0,
                            "Increased in AD",
                            "Decreased in AD"))

plot_df_uncorrected <- plot_df_uncorrected %>%
  left_join(direction_map, by = "Genus") %>%
  mutate(RA_adj = RA + 1e-6)

longgenus$RA_adj <- longgenus$RA + 1e-6
prevalence <- longgenus %>%
  group_by(Genus) %>%
  summarize(prev = sum(RA > 0) / n())
keep_genera <- prevalence %>% filter(prev >= 0.2) %>% pull(Genus)
mean_abundance <- longgenus %>%
  group_by(Genus) %>%
  summarize(meanRA = mean(RA_adj))
keep_abundant <- mean_abundance %>% filter(meanRA > 0.001) %>% pull(Genus)
filtered_genera <- intersect(keep_genera, keep_abundant)
filtered_df <- plot_df_uncorrected %>% 
  filter(Genus %in% filtered_genera)
plot_inc_filt <- filtered_df %>% filter(direction == "Increased in AD")
plot_dec_filt <- filtered_df %>% filter(direction == "Decreased in AD")
plot_inc_filt$diagnosis_clean <- recode(plot_inc_filt$diagnosis,
                                        "stool sample AD" = "AD",
                                        "stool sample Control" = "Control")

plot_dec_filt$diagnosis_clean <- recode(plot_dec_filt$diagnosis,
                                        "stool sample AD" = "AD",
                                        "stool sample Control" = "Control")
ggplot(plot_inc_filt,
       aes(x = diagnosis_clean, y = RA_adj, fill = Phylum)) +
  geom_boxplot(alpha = 0.9, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
  scale_y_log10() +
  annotation_logticks(sides = "l") +
  facet_wrap(~ Genus, ncol = 4) +  
  theme_bw(base_size = 14) +
  labs(
    title = "Genera Increased in AD (Uncorrected Wilcoxon p < 0.05)",
    y = "Relative Abundance (log10)",
    x = ""
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "right",
    panel.spacing = unit(0.8, "lines"),
    plot.title = element_text(face = "bold")
  )


plot_df_uncorrected$Genus_wrapped <- str_wrap(plot_df_uncorrected$Genus, width = 10)
plot_inc <- plot_df_uncorrected %>% filter(direction == "Increased in AD")
plot_dec <- plot_df_uncorrected %>% filter(direction == "Decreased in AD")

plot_inc$diagnosis_clean <- recode(plot_inc$diagnosis,
                                   "stool sample AD" = "AD",
                                   "stool sample Control" = "Control")

plot_dec$diagnosis_clean <- recode(plot_dec$diagnosis,
                                   "stool sample AD" = "AD",
                                   "stool sample Control" = "Control")
top15 <- longgenus %>%
  group_by(Genus) %>%
  summarise(total_abundance = sum(RA, na.rm = TRUE)) %>%
  arrange(desc(total_abundance)) %>%
  slice(1:15) %>%
  pull(Genus)

plot_inc_top15 <- plot_inc %>% filter(Genus %in% top15)
plot_dec_top15 <- plot_dec %>% filter(Genus %in% top15)


ggplot(plot_inc_top15, 
       aes(x = diagnosis_clean, y = RA_adj, fill = Phylum)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 1) +
  scale_y_log10() +
  annotation_logticks(sides = "l") +
  facet_wrap(~ Genus, ncol = 4) +
  theme_bw(base_size = 14) +
  labs(
    title = "Out of the Top 15, Genera More Abundant in AD (Uncorrected p < 0.05)",
    y = "Relative Abundance (log10 scale)",
    x = ""
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "right",
    panel.spacing = unit(0.8, "lines"),
    plot.title = element_text(face = "bold")
  )

ggplot(plot_dec_top15, 
       aes(x = diagnosis_clean, y = RA_adj, fill = Phylum)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 1) +
  scale_y_log10() +
  annotation_logticks(sides = "l") +
  facet_wrap(~ Genus, ncol = 4) +
  theme_bw(base_size = 14) +
  labs(
    title = "Out of the Top 15, Genera Less Abundant in AD (Uncorrected p < 0.05)",
    y = "Relative Abundance (log10 scale)",
    x = ""
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "right",
    panel.spacing = unit(0.8, "lines"),
    plot.title = element_text(face = "bold")
  )
#### differential abundance analysis
#####filtering and factoring

meta_ps$dx <- case_when(
  meta_ps$diagnosis == "stool sample AD" ~ "AD",
  meta_ps$diagnosis == "stool sample Control" ~ "Control",
  TRUE ~ meta_ps$diagnosis
)

meta_ps$dx <- factor(meta_ps$dx, levels = c("Control", "AD"))

longgenus$RA_adj <- longgenus$RA + 1e-6

prev_df <- longgenus %>%
  group_by(Genus) %>%
  summarise(prev = sum(RA > 0) / n())

abund_df <- longgenus %>%
  group_by(Genus) %>%
  summarise(meanRA = mean(RA_adj))

filtered_genera <- prev_df %>%
  filter(prev >= 0.2) %>%
  inner_join(abund_df, by = "Genus") %>%
  filter(meanRA > 0.001) %>%
  pull(Genus)

genus_mat_full <- asv.genus.sum %>%
  column_to_rownames("Genus") %>%
  as.matrix()

genus_mat_filt <- genus_mat_full[rownames(genus_mat_full) %in% filtered_genera, ]
genus_mat_filt <- genus_mat_filt[, meta_ps$sample_id]  # match metadata



ps_genus <- phyloseq(
  otu_table(genus_mat_filt, taxa_are_rows = TRUE),
  sample_data(meta_ps)
)
#####correct diiferential abundance
dds <- phyloseq_to_deseq2(ps_genus, ~ dx)
dds <- DESeq(dds, fitType = "local")
res <- results(dds, contrast = c("dx", "AD", "Control"))
res_df <- as.data.frame(res) %>%
  rownames_to_column("Genus") %>%
  arrange(padj) %>%
  mutate(
    neg_log10_padj = -log10(padj),
    sig = ifelse(!is.na(padj) & padj < 0.05, "sig", "ns"),
    direction = case_when(
      padj < 0.05 & log2FoldChange > 0 ~ "Increased in AD",
      padj < 0.05 & log2FoldChange < 0 ~ "Decreased in AD",
      TRUE ~ "NS"
    )
  )

sig_genera <- res_df %>% filter(padj < 0.05) %>% pull(Genus)

#### volcano
ggplot(res_df, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(aes(color = direction), size = 3, alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  scale_color_manual(
    values = c(
      "Increased in AD"   = "firebrick",
      "Decreased in AD"   = "steelblue",
      "NS"                = "grey60"
    )
  ) +
  theme_bw(base_size = 14) +
  labs(
    title = "Volcano Plot showing AD vs Control (Filtered DESeq2)",
    x = "log2 Fold Change (AD / Control)",
    y = "-log10 Adjusted p-value",
    color = ""
  ) +
  theme(
    plot.title = element_text(face = "bold")
  )

#####boxplot significant filetered
plot_df_sig <- longgenus %>%
  filter(Genus %in% sig_genera) %>%
  mutate(
    diagnosis_clean = recode(diagnosis,
                             "stool sample AD" = "AD",
                             "stool sample Control" = "Control"),
    RA_adj = RA + 1e-6
  )

plotsig<- ggplot(plot_df_sig,
       aes(x = diagnosis_clean, y = RA_adj, fill = diagnosis_clean)) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
  facet_wrap(~ Genus, scales = "free_y") +
  scale_y_log10() +
  theme_bw(base_size = 14) +
  labs(
    title = "DESeq2 Significantly Differential Genera (Filtered)",
    y = "Relative Abundance (log10)",
    x = ""
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )
ggsave("deseq.genera.2.png", plot = plotsig, width = 14, height = 6, units = "in")

####without filtering
genus_mat_unfilt <- asv.genus.sum %>%
  column_to_rownames("Genus") %>%
  as.matrix()
genus_mat_unfilt <- genus_mat_unfilt[, meta_ps$sample_id]

ps_genus_unfilt <- phyloseq(
  otu_table(genus_mat_unfilt, taxa_are_rows = TRUE),
  sample_data(meta_ps)
)
dds_unfilt <- phyloseq_to_deseq2(ps_genus_unfilt, ~ dx)
dds_unfilt <- DESeq(dds_unfilt, fitType = "local")
res_unfilt <- results(dds_unfilt, contrast = c("dx", "AD", "Control"))
res_df_unfilt <- as.data.frame(res_unfilt) %>%
  rownames_to_column("Genus")
sig_genera_unfilt <- res_df_unfilt %>%
  filter(padj < 0.05) %>%
  pull(Genus)

plot_df_sig_unfilt <- longgenus %>%
  filter(Genus %in% sig_genera_unfilt) %>%
  mutate(
    diagnosis_clean = recode(
      diagnosis,
      "stool sample AD" = "AD",
      "stool sample Control" = "Control"
    ),
    RA_adj = RA + 1e-6
  )

plotsig_unfilt <- ggplot(plot_df_sig_unfilt,
                         aes(x = diagnosis_clean, y = RA_adj, fill = diagnosis_clean)) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
  facet_wrap(~ Genus, scales = "free_y") +
  scale_y_log10() +
  theme_bw(base_size = 14) +
  labs(
    title = "DESeq2 Significantly Differential Genera (Unfiltered)",
    y = "Relative Abundance (log10)",
    x = ""
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )
ggsave("deseq.genera.unfiltered.png",
       plot = plotsig_unfilt,
       width = 16, height = 8, units = "in")

#barplot 
res_df_sorted <- res_df %>% arrange(log2FoldChange)

ggplot(res_df_sorted,
       aes(x = reorder(Genus, log2FoldChange),
           y = log2FoldChange,
           fill = log2FoldChange > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c("TRUE" = "firebrick", "FALSE" = "steelblue")
  ) +
  theme_bw(base_size = 14) +
  labs(
    title = "Log2 Fold Change (Filtered Genera)",
    x = "Genus",
    y = "log2 Fold Change (AD / Control)",
    fill = "Higher in AD?"
  ) +
  theme(plot.title = element_text(face = "bold"))


##uncorrected significant differeces using wilcoxon
genus_tests_uncorrected <- longgenus %>%
  group_by(Genus) %>%
  wilcox_test(RA ~ diagnosis) %>%
  select(Genus, p)

sig_uncorrected_genera <- genus_tests_uncorrected %>%
  filter(p < 0.05) %>%
  pull(Genus)

plot_df_uncorrected_sig <- longgenus %>%
  filter(Genus %in% sig_uncorrected_genera) %>%
  mutate(
    diagnosis_clean = recode(
      diagnosis,
      "stool sample AD" = "AD",
      "stool sample Control" = "Control"
    ),
    RA_adj = RA + 1e-6    
  )
plot_uncorrected <- ggplot(plot_df_uncorrected_sig,
                           aes(x = diagnosis_clean, y = RA_adj, fill = diagnosis_clean)) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
  facet_wrap(~ Genus, scales = "free_y") +
  scale_y_log10() +
  theme_bw(base_size = 14) +
  labs(
    title = "Uncorrected Wilcoxon Significant Genera (p < 0.05)",
    y = "Relative Abundance (log10)",
    x = ""
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )
ggsave("uncorrected.wilcoxon.genera.png",
       plot = plot_uncorrected,
       width = 24, height = 12, units = "in")

##uncorrected significant differeces using wilcoxon with filtering 
####preprocessing
longgenus$RA_adj <- longgenus$RA + 1e-6
prev_df_filt <- longgenus %>%
  group_by(Genus) %>%
  summarise(prev = sum(RA > 0) / n())
abund_df_filt <- longgenus %>%
  group_by(Genus) %>%
  summarise(meanRA = mean(RA_adj))

filtered_genera <- prev_df_filt %>%
  filter(prev >= 0.2) %>%
  inner_join(abund_df_filt, by = "Genus") %>%
  filter(meanRA > 0.001) %>%
  pull(Genus)

genus_tests_filtered <- longgenus %>%
  filter(Genus %in% filtered_genera) %>%
  group_by(Genus) %>%
  wilcox_test(RA ~ diagnosis) %>%
  select(Genus, p)

sig_uncorrected_filt <- genus_tests_filtered %>%
  filter(p < 0.05) %>%
  pull(Genus)
plot_df_uncorrected_filt <- longgenus %>%
  filter(Genus %in% sig_uncorrected_filt) %>%
  mutate(
    diagnosis_clean = recode(
      diagnosis,
      "stool sample AD" = "AD",
      "stool sample Control" = "Control"
    ),
    RA_adj = RA + 1e-6
  )

plot_uncorrected_filtered <- ggplot(plot_df_uncorrected_filt,
                                    aes(x = diagnosis_clean, y = RA_adj, fill = diagnosis_clean)) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
  facet_wrap(~ Genus, scales = "free_y") +
  scale_y_log10() +
  theme_bw(base_size = 14) +
  labs(
    title = "Uncorrected Wilcoxon Significant Genera (Filtered Dataset)",
    y = "Relative Abundance (log10 scale)",
    x = ""
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )
ggsave("uncorrected.wilcoxon.filtered.genera.png",
       plot = plot_uncorrected_filtered,
       width = 26, height = 14, units = "in")

######### new filtered corrected boxplots but with less in AD vs. more in AD
sig_res <- res_df %>% 
  filter(Genus %in% sig_genera) %>%
  mutate(direction = ifelse(log2FoldChange > 0, "Increased in AD", "Decreased in AD"))

sig_inc <- sig_res %>% filter(direction == "Increased in AD") %>% pull(Genus)
sig_dec <- sig_res %>% filter(direction == "Decreased in AD") %>% pull(Genus)

plot_df_sig_filt <- longgenus %>%
  mutate(
    diagnosis_clean = recode(
      diagnosis,
      "stool sample AD" = "AD",
      "stool sample Control" = "Control"
    ),
    RA_adj = RA + 1e-6
  )

plot_inc <- plot_df_sig_filt %>%
  filter(Genus %in% sig_inc)

plot_inc_fig <- ggplot(plot_inc,
                       aes(x = diagnosis_clean, y = RA_adj, fill = diagnosis_clean)) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
  facet_wrap(~ Genus, scales = "free_y") +
  scale_y_log10() +
  theme_bw(base_size = 14) +
  labs(
    title = "DESeq2 Significant Genera Increased in AD (Filtered)",
    y = "Relative Abundance (log10)",
    x = ""
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave("deseq.sig.increased.AD.filtered.png",
       plot = plot_inc_fig,
       width = 14, height = 10, units = "in")

plot_dec <- plot_df_sig_filt %>%
  filter(Genus %in% sig_dec)

plot_dec_fig <- ggplot(plot_dec,
                       aes(x = diagnosis_clean, y = RA_adj, fill = diagnosis_clean)) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
  facet_wrap(~ Genus, scales = "free_y") +
  scale_y_log10() +
  theme_bw(base_size = 14) +
  labs(
    title = "DESeq2 Significant Genera Decreased in AD (Filtered)",
    y = "Relative Abundance (log10)",
    x = ""
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave("deseq.sig.decreased.AD.filtered.png",
       plot = plot_dec_fig,
       width = 18, height = 10, units = "in")
