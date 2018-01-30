### campy w/ variable diet - final microbiome analyses

library(ggplot2)
library(phyloseq)
library(reshape2)
library(vegan)
SEED_INT = 10
set.seed(SEED_INT) # set the seed number for reproducibility in randomly generated numbers

### DATA INPUT AND RESHAPING ------------------------------------------------
# read in the output file from DADA2 and metadata
ps <- readRDS("~/Documents/projects/campy_murine_diets/data/campy_phyloseq_obj_less_strict.rds")
# remove neg control
neg_control = prune_samples(rownames(sample_data(ps)) %in% c("Plate3-D03"),ps)
path = "~/Documents/projects/campy_murine_diets/data/"
metadata = read.table(paste0(path,"Campy3microbiome1-23-2017.txt"),header=TRUE,
                      sep = "\t")
colnames(metadata)[colnames(metadata) == 'Plate'] = "plate_sample"
# read in mouse weights and add to the metadata
weights = read.table(paste0(path,"campystudy_weights.txt"),header=TRUE,
                     sep = "\t",check.names=FALSE)
weights_long = melt(weights,id.vars=c('Day'))
colnames(weights_long)[colnames(weights_long) == "variable"] = "Sample"
colnames(weights_long)[colnames(weights_long) == "value"] = "weight"
colnames(weights_long)[colnames(weights_long) == "Day"] = "days_post_infection"
metadata = merge(metadata,weights_long,all=TRUE)
#metadata = merge(metadata,shedding_long,all=TRUE)
# filter out samples with plate_sample as NA -- these don't have microbiome sequencing
metadata = metadata[!is.na(metadata$plate_sample),]

# merge metadata with phyloseq sample data
ps_sampledata_asdf = data.frame(sample_data(ps))
# There is a sample with sequencing data thhat is not in the sample sheet. Filter out for now.
rownames(metadata) = metadata$plate_sample
ps_sampledata_asdf = ps_sampledata_asdf[rownames(ps_sampledata_asdf) %in% rownames(metadata),]

samples = rownames(ps_sampledata_asdf)
# filter metadata by samples that have sequencing data
filtered_metadata = metadata[rownames(metadata) %in% samples,]
new_sampledata = merge(filtered_metadata,ps_sampledata_asdf,by="row.names")
rownames(new_sampledata) = new_sampledata$Row.names
new_sampledata$Row.names = NULL
# subset the samples in the phyloseq object, then reassign the metadata as sample_data
ps = prune_samples(rownames(sample_data(ps)) %in% rownames(new_sampledata), ps)
#sample_data(ps) = new_sampledata
sample_data(ps) = new_sampledata

# For this study, remove dN samples
ps = prune_samples(sample_data(ps)$diet %in% c("HC","dPD","dZD"),ps)
# remove samples that have 0 reads
ps = prune_samples(sample_sums(ps) > 0, ps)

###--------------------------------------------------------------------------

### OTU AND SAMPLE TRIMMING--------------------------------------------------
# prune samples with fewer than 10 counts, and those with fewer than 5 non-zero OTUS
old_samples = rownames(sample_data(ps))
ps_trimmed = prune_samples(sample_sums(ps)>=10, ps)
# Remove OTUs with 20 or fewer counts summed across all samples
ps_trimmed = prune_taxa(colSums(otu_table(ps_trimmed)) > 20,ps_trimmed)
# get samples with fewer than 5 non-zero OTUs
nonzero_count = apply(otu_table(ps_trimmed), 1, function(c)sum(c!= 0))

ps_trimmed = prune_samples(nonzero_count > 5, ps_trimmed)



# keep track of the trimmed samples and taxa
pruned_samples = old_samples[!old_samples %in% rownames(sample_data(ps_trimmed))]
old_taxa = rownames(tax_table(ps))

# generate a phyloseq object converted to relative abundance instead of counts
ps_trimmed_relative  = transform_sample_counts(ps_trimmed, function(x) x / sum(x) )


###--------------------------------------------------------------------------

### FIGURE 1: alpha diversity over time for each diet and each infection status
# We will rarefy to 5000 counts per sample, which was the approximate point at which richness plateaued
# when visualized with ggrare

# for alpha diversity, we need to used the untrimmed otu table, with the exception of 
# removing the samples that had 0 reads.
rare_min = 5000
# remove samples with fewer counts than the minimum
ps_rare = prune_samples(sample_sums(ps) >=rare_min,ps)
ps_rare = rarefy_even_depth(ps_rare,sample.size=5000,rngseed=SEED_INT) # OTUs no longer in the dataset will be removed


days_to_plot = c(0,1,5,9) # select which days we'd like to plot and compare
ps_for_alpha = prune_samples(sample_data(ps_rare)$days_post_infection %in% days_to_plot,ps_rare)
# prune the taxa that don't appear in any of these samples because they were only present in excluded time points
ps_for_alpha <- filter_taxa(ps_for_alpha, function(x) sum(x > 1) > 0, TRUE)

sample_data(ps_for_alpha)$days_post_infection = as.factor(sample_data(ps_for_alpha)$days_post_infection)
# get actual values and test differences, evenness first
evenness = estimate_richness(ps_for_alpha,measures="InvSimpson")
asdf = as(sample_data(ps_for_alpha),"data.frame")
asdf = merge(asdf,evenness,by="row.names")
asdf$days_post_infection = asdf$days_post_infection
asdf$infected = as.factor(asdf$infected)
diet_anova = aov(InvSimpson~diet*days_post_infection*infected,data=asdf)
diet_tukey_posthoc_evenness = TukeyHSD(diet_anova)

# repeat for richness
richness = estimate_richness(ps_for_alpha,measures="Observed")
# for some reason row names get mangled when using measures="Observed")...
rownames(richness) = gsub("\\.","-",rownames(richness))
asdf = as(sample_data(ps_for_alpha),"data.frame")
asdf = merge(asdf,richness,by="row.names")
asdf$days_post_infection = asdf$days_post_infection
asdf$infected = as.factor(asdf$infected)
diet_anova = aov(Observed~diet*days_post_infection*infected,data=asdf)
diet_tukey_posthoc_richness = TukeyHSD(diet_anova)

# plot evenness and save associated test object
theme_set(theme_bw())
p = plot_richness(ps_for_alpha,measures = "Observed", x='days_post_infection')
p = p + geom_boxplot(alpha=0.5) + facet_grid(infected~diet)
#p$layers = p$layers[-1]
p
# save as svg
ggsave(file="~/Documents/projects/campy_murine_diets/results/Observed_div_all.svg", plot=p, width=10, height=8)
# save the anova results
saveRDS(diet_tukey_posthoc_evenness, "~/Documents/projects/campy_murine_diets/results/Observed_anova_posthoc.rds")

# plot richness and save associated test object
theme_set(theme_bw())
p = plot_richness(ps_for_alpha,measures = "InvSimpson", x='days_post_infection')
p = p + geom_boxplot(alpha=0.5) + facet_grid(infected~diet)
#p$layers = p$layers[-1]
p
# save as svg
ggsave(file="~/Documents/projects/campy_murine_diets/results/InvSimpson_div_all.svg", plot=p, width=10, height=8)
# save the anova results
saveRDS(diet_tukey_posthoc_evenness, "~/Documents/projects/campy_murine_diets/results/InvSimpson_anova_posthoc.rds")


###------------------------------------------------------------------------

### correlation between campy shedding and alpha diversity

copy_sample_data = sample_data(ps_for_alpha)
richness$plate_sample = rownames(richness)
evenness$plate_sample = rownames(evenness)

copy_sample_data = merge(copy_sample_data,richness,by='plate_sample')

# plot richness vs. campy for all samples
copy_sample_data

# plot evenness vs. campy for all samples



###------------------------------------------------------------------------

### FIGURE 2: NMDS of beta diversity for all samples-----------------------
days_to_plot = c(0,1,5,9)
ps_for_beta = prune_samples(sample_data(ps_trimmed)$days_post_infection %in% days_to_plot,ps_trimmed)
sample_data(ps_for_beta)$days_post_infection = as.factor(sample_data(ps_for_beta)$days_post_infection)

ps_trimmed.ord <- ordinate(ps_for_beta, "NMDS", "bray")
theme_set(theme_bw())
p1 = plot_ordination(ps_for_beta, ps_trimmed.ord, color="days_post_infection",type="samples", shape="infected") + facet_wrap(~diet) + theme(strip.text = element_text(size=20)) 
p1 = p1 + geom_point(size=5,alpha=0.5)
#new_data = newdata=data.frame(ps_d015.ord$points[,1],ps_d015.ord$points[,2],sample_data(ps_d015)$infected)
#colnames(new_data) = c('NMDS1','NMDS2','infected')
#p1 = p1 + geom_point(data=new_data,mapping=aes(shape='infected'))
p1
ggsave(file="~/Documents/projects/campy_murine_diets/results/beta_div_bray_all.svg", plot=p1, width=15, height=4.5)
###------------------------------------------------------------------------

### NMDS within each diet--------------------------------------------------
days_to_plot = c(0,1,5,9)
ps_for_beta = prune_samples(sample_data(ps_trimmed)$days_post_infection %in% days_to_plot,ps_trimmed)
sample_data(ps_for_beta)$days_post_infection = as.factor(sample_data(ps_for_beta)$days_post_infection)

# subset by diet
ps_for_beta_dPD = prune_samples(sample_data(ps_for_beta)$diet == "dPD",ps_for_beta)
ps_for_beta_dZD = prune_samples(sample_data(ps_for_beta)$diet == "dZD",ps_for_beta)
ps_for_beta_HC = prune_samples(sample_data(ps_for_beta)$diet == "HC",ps_for_beta)

ps_trimmed_dPD.ord <- ordinate(ps_for_beta_dPD, "NMDS", "bray")
ps_trimmed_dZD.ord <- ordinate(ps_for_beta_dZD, "NMDS", "bray")
ps_trimmed_HC.ord <- ordinate(ps_for_beta_HC, "NMDS", "bray")

# perform adonis test within each diet to 
phyloseq::distance(erie_scale, method = "bray")
adonis(ps_trimmed_dPD.ord,)

theme_set(theme_bw())
p1 = plot_ordination(ps_for_beta_dPD, ps_trimmed_dPD.ord, color="days_post_infection",type="samples", shape="infected") + theme(strip.text = element_text(size=20)) 
p1 = p1 + geom_point(size=5,alpha=0.5)
#new_data = newdata=data.frame(ps_d015.ord$points[,1],ps_d015.ord$points[,2],sample_data(ps_d015)$infected)
#colnames(new_data) = c('NMDS1','NMDS2','infected')
#p1 = p1 + geom_point(data=new_data,mapping=aes(shape='infected'))
p1
ggsave(file="~/Documents/projects/campy_murine_diets/results/beta_div_bray_dPD.svg", plot=p1, width=5, height=3)


p1 = plot_ordination(ps_for_beta_dZD, ps_trimmed_dZD.ord, color="days_post_infection",type="samples", shape="infected") + theme(strip.text = element_text(size=20)) 
p1 = p1 + geom_point(size=5,alpha=0.5)
#new_data = newdata=data.frame(ps_d015.ord$points[,1],ps_d015.ord$points[,2],sample_data(ps_d015)$infected)
#colnames(new_data) = c('NMDS1','NMDS2','infected')
#p1 = p1 + geom_point(data=new_data,mapping=aes(shape='infected'))
p1
ggsave(file="~/Documents/projects/campy_murine_diets/results/beta_div_bray_dZD.svg", plot=p1, width=5, height=3)


p1 = plot_ordination(ps_for_beta_HC, ps_trimmed_HC.ord, color="days_post_infection",type="samples", shape="infected") + theme(strip.text = element_text(size=20)) 
p1 = p1 + geom_point(size=5,alpha=0.5)
#new_data = newdata=data.frame(ps_d015.ord$points[,1],ps_d015.ord$points[,2],sample_data(ps_d015)$infected)
#colnames(new_data) = c('NMDS1','NMDS2','infected')
#p1 = p1 + geom_point(data=new_data,mapping=aes(shape='infected'))
p1
ggsave(file="~/Documents/projects/campy_murine_diets/results/beta_div_bray_HC.svg", plot=p1, width=5, height=3)

###------------------------------------------------------------------------


### DIFFERENTIAL ABUNDANCE-------------------------------------------------
library("DESeq2")
alpha = 0.05

### Compare each diet at d0 first, since microbiota might be implicated in susceptibility
ps_d0 = prune_samples(sample_data(ps_trimmed)$days_post_infection %in% c(0),ps_trimmed)
# trim again to remove low-abundance taxa within this subset. Loosen the # samples restriction
# since this is a smaller sample size
ps_d0_taxatrim = filter_taxa(ps_d0, function(x) sum(x > 1) > 3, TRUE)
ps_d0_taxatrim = filter_taxa(ps_d0_taxatrim, function(x) max(x) > 10, TRUE)

# subset by diet to make pair-wise comparisons
pair_diet_comparison = function(x,y,path,alpha) {
  comp = prune_samples(sample_data(ps_d0_taxatrim)$diet %in% c(x,y),ps_d0_taxatrim)
  comp = filter_taxa(comp, function(x) sum(x > 1) > 3, TRUE)
  comp_deseq = phyloseq_to_deseq2(comp, ~ diet)
  comp_deseq = DESeq(comp_deseq, test="Wald", fitType="parametric")
  comp_deseq = results(comp_deseq)
  # replace NA padj w/ 0
  comp_deseq$padj[is.na(comp_deseq$padj)] <- 0
  comp_deseq = comp_deseq[comp_deseq$padj < alpha,]
  comp_deseq = merge(comp_deseq,tax_table(comp)[rownames(comp_deseq),],by="row.names")
  # save
  write.table(comp_deseq, file = path, sep = "\t", quote = FALSE,row.names = FALSE)
  return(comp_deseq)
}
#d0_dNvdPD = pair_diet_comparison(x="dN",y="dPD",path="~/Documents/projects/campy_murine_diets/results/d0_dNvdPD_diff.txt",alpha=0.05)
#d0_dNvdZD = pair_diet_comparison(x="dN",y="dZD",path="~/Documents/projects/campy_murine_diets/results/d0_dNvdZD_diff.txt",alpha=0.05)
#d0_dNvHC = pair_diet_comparison(x="dN",y="HC",path="~/Documents/projects/campy_murine_diets/results/d0_dNvHC_diff.txt",alpha=0.05)
d0_dZDvdPD = pair_diet_comparison(x="dZD",y="dPD",path="~/Documents/projects/campy_murine_diets/results/d0_dZDvdPD_diff.txt",alpha=0.05)
d0_dZDvHC = pair_diet_comparison(x="dZD",y="HC",path="~/Documents/projects/campy_murine_diets/results/d0_dZDvHC_diff.txt",alpha=0.05)
d0_dPDvHC = pair_diet_comparison(x="dPD",y="HC",path="~/Documents/projects/campy_murine_diets/results/d0_dPDvHC_diff.txt",alpha=0.05)

# difference between dZD and the rest?
higher_in_dZD = d0_dZDvdPD[d0_dZDvdPD$log2FoldChange > 0,]
higher_in_dZD = d0_dZDvHC[d0_dZDvHC$log2FoldChange > 0,]
#higher_in_dZD = d0_dNvdZD[d0_dNvdZD$log2FoldChange < 0,]

lower_in_dZD = d0_dZDvdPD[d0_dZDvdPD$log2FoldChange < 0,][c("Row.names","Order","Family","Genus")]
lower_in_dZD = merge(lower_in_dZD,d0_dZDvHC[d0_dZDvHC$log2FoldChange < 0,][c("Row.names","Order","Family","Genus")],by="Row.names")
#lower_in_dZD = merge(lower_in_dZD,d0_dNvdZD[d0_dNvdZD$log2FoldChange > 0,][c("Row.names","Order","Family","Genus")],by="Row.names")

### Now compare infected vs. uninfected at each time point for each diet
# Effect of infection at each time point, for each diet
ps_d1 = prune_samples(sample_data(ps_trimmed)$days_post_infection %in% c(1),ps_trimmed)
ps_d1_taxatrim = filter_taxa(ps_d1, function(x) sum(x > 1) > 3, TRUE)
ps_d1_taxatrim = filter_taxa(ps_d1_taxatrim, function(x) max(x) > 10, TRUE)
ps_d5 = prune_samples(sample_data(ps_trimmed)$days_post_infection %in% c(5),ps_trimmed)
ps_d5_taxatrim = filter_taxa(ps_d5, function(x) sum(x > 1) > 3, TRUE)
ps_d5_taxatrim = filter_taxa(ps_d5_taxatrim, function(x) max(x) > 10, TRUE)
ps_d9 = prune_samples(sample_data(ps_trimmed)$days_post_infection %in% c(9),ps_trimmed)
ps_d9_taxatrim = filter_taxa(ps_d9, function(x) sum(x > 1) > 3, TRUE)
ps_d9_taxatrim = filter_taxa(ps_d9_taxatrim, function(x) max(x) > 10, TRUE)

infection_comparison = function(x,ps_obj,path,alpha) {
  comp = prune_samples(sample_data(ps_obj)$diet %in% c(x),ps_obj)
  comp = filter_taxa(comp, function(x) sum(x > 1) > 3, TRUE)
  
  comp_deseq = phyloseq_to_deseq2(comp, ~ infected)
  print(comp_deseq)
  comp_deseq = DESeq(comp_deseq, test="Wald", fitType="parametric")
  comp_deseq = results(comp_deseq)
  print(sum(is.na(comp_deseq)))
  comp_deseq = merge(comp_deseq,tax_table(comp)[rownames(comp_deseq),],by="row.names")
  print(sum(is.na(comp_deseq)))
  comp_deseq = comp_deseq[comp_deseq$padj < alpha,]
  # save
  write.table(comp_deseq, file = path, sep = "\t", quote = FALSE,row.names = FALSE)
  return(comp_deseq)
}
# run the comparison for all time points and save the results.
#d1_dN = infection_comparison(x="dN",ps_obj=ps_d1_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d1_dN_infection_diff.txt",alpha=0.05)
#d5_dN = infection_comparison(x="dN",ps_obj=ps_d5_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d5_dN_infection_diff.txt",alpha=0.05)
#d9_dN = infection_comparison(x="dN",ps_obj=ps_d9_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d9_dN_infection_diff.txt",alpha=0.05)
d1_dPD = infection_comparison(x="dPD",ps_obj=ps_d1_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d1_dPD_infection_diff.txt",alpha=0.05)
d5_dPD = infection_comparison(x="dPD",ps_obj=ps_d5_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d5_dPD_infection_diff.txt",alpha=0.05)
d9_dPD = infection_comparison(x="dPD",ps_obj=ps_d9_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d9_dPD_infection_diff.txt",alpha=0.05)
d1_dZD = infection_comparison(x="dZD",ps_obj=ps_d1_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d1_dZD_infection_diff.txt",alpha=0.05)
d5_dZD = infection_comparison(x="dZD",ps_obj=ps_d5_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d5_dZD_infection_diff.txt",alpha=0.05)
d9_dZD = infection_comparison(x="dZD",ps_obj=ps_d9_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d9_dZD_infection_diff.txt",alpha=0.05)
d1_HC = infection_comparison(x="HC",ps_obj=ps_d1_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d1_HC_infection_diff.txt",alpha=0.05)
d5_HC = infection_comparison(x="HC",ps_obj=ps_d5_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d5_HC_infection_diff.txt",alpha=0.05)
d9_HC = infection_comparison(x="HC",ps_obj=ps_d9_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d9_HC_infection_diff.txt",alpha=0.05)

# at each time point, determine the shared differentially abundant sequences for each diet
# day 1, up in all diets all diets
d1_dN[d1_dN$log2FoldChange > 0,]

###------------------------------------------------------------------------




