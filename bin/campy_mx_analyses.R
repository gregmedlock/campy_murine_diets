### campy w/ variable diet - final microbiome analyses

library(ggplot2)
library(phyloseq)
library(reshape2)
library(vegan)
library(gridExtra)
SEED_INT = 10
set.seed(SEED_INT) # set the seed number for reproducibility in randomly generated numbers

### DATA INPUT AND RESHAPING ------------------------------------------------
# read in the output file from DADA2 and metadata
ps <- readRDS("~/Documents/projects/campy_murine_diets/data/campy_phyloseq_obj_less_strict.rds")
# remove neg control
neg_control = prune_samples(rownames(sample_data(ps)) %in% c("Plate3-D03"),ps)
# save a bar plot of the neg control (mock) for quality assurance supplementary figure
p = plot_bar(neg_control, "Order", fill="Order")
p
ggsave(file="~/Documents/projects/campy_murine_diets/results/mock_abundances.svg",p,height=10,width=15)

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

# remove mitochondrial and chloroplast reads
og_taxa = data.frame(tax_table(ps))
# first remove non-bacteria
keep = rownames(og_taxa[og_taxa$Kingdom == "Bacteria",])
ps = prune_taxa(keep,ps)
keep = rownames(og_taxa[og_taxa$Family != "Mitochondria",])
ps = prune_taxa(keep,ps)
keep = rownames(og_taxa[og_taxa$Class != "Chloroplast",])
ps = prune_taxa(keep,ps)

# For this study, remove dN samples
ps = prune_samples(sample_data(ps)$diet %in% c("HC","dPD","dZD"),ps)
# remove samples that have 0 reads
ps = prune_samples(sample_sums(ps) > 0, ps)
# get rid of day 11
ps = prune_samples(sample_data(ps)$days_post_infection %in% c(0,1,5,9),ps)


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

# define variables for ggplot sizes to make them consistent across figures
legend_textsize=20
legend_symbolsize=1
size_striptext = 20
size_title = 20
size_axistext = 12

# plot evenness and save associated test object
theme_set(theme_bw())
infection_labeller = list("FALSE"="Uninfected","TRUE"="Infected")
# change infected label from True/False to infected/uninfected
tempsampledata = sample_data(ps_for_alpha)
col2=ifelse(tempsampledata$infected==TRUE,"Infected","Uninfected")
tempsampledata$infected = col2
sample_data(ps_for_alpha) = tempsampledata
p = plot_richness(ps_for_alpha,measures = "Observed", x='days_post_infection',color='days_post_infection')
p1 = p + geom_boxplot(alpha=0.5) + facet_grid(infected~diet) +
  theme(strip.text.x=element_text(size=size_striptext,face="bold"),
        strip.text.y=element_text(size=size_striptext,face="bold"),
        axis.title.x = element_text(size=size_title), 
        axis.text.x = element_text(size=size_axistext,angle=0,hjust=0.5),
        axis.text.y = element_text(size=size_axistext),
        axis.title.y = element_text(size=size_title),
        legend.position="none",
        plot.margin = unit(c(1,1,1,1), "cm")) + 
        labs(y="Richness (observed ASVs)",x="days post-infection")
        
#p$layers = p$layers[-1]
p1
# save as svg
ggsave(file="~/Documents/projects/campy_murine_diets/results/Observed_div_all.svg", plot=p1, width=10, height=6)
# save the anova results
saveRDS(diet_tukey_posthoc_richness, "~/Documents/projects/campy_murine_diets/results/Observed_anova_posthoc.rds")

# plot richness and save associated test object
theme_set(theme_bw())
p = plot_richness(ps_for_alpha,measures = "InvSimpson", x='days_post_infection', color='days_post_infection')
p2 = p + geom_boxplot(alpha=0.5) + facet_grid(infected~diet) +
  theme(strip.text.x=element_text(size=size_striptext,face="bold"),
        strip.text.y=element_text(size=size_striptext,face="bold"),
        axis.title.x = element_text(size=size_title), 
        axis.text.x = element_text(size=size_axistext,angle=0,hjust=0.5),
        axis.text.y = element_text(size=size_axistext),
        axis.title.y = element_text(size=size_title),
        legend.position="none",
        plot.margin = unit(c(1,1,1,1), "cm")) + 
        labs(y="Evenness (Inverse Simpson)",x="days post-infection") +
        guides(colour = guide_legend(override.aes = list(size=legend_symbolsize)))
#p$layers = p$layers[-1]
p2
# save as svg
ggsave(file="~/Documents/projects/campy_murine_diets/results/InvSimpson_div_all.svg", plot=p2, width=10, height=6)
# save the anova results
saveRDS(diet_tukey_posthoc_evenness, "~/Documents/projects/campy_murine_diets/results/InvSimpson_anova_posthoc.rds")

###------------------------------------------------------------------------

### correlation between campy shedding and alpha diversity

copy_sample_data = sample_data(ps_for_alpha)
richness$plate_sample = rownames(richness)
evenness$plate_sample = rownames(evenness)

copy_sample_data = merge(data.frame(copy_sample_data),richness,by='plate_sample')
copy_sample_data = merge(copy_sample_data,evenness,by='plate_sample')

# generate expression to italicize campy in labels
xlabel_expression = expression(paste(italic("C. jejuni"),"/10mg stool"))
# plot richness vs. campy for all samples
# plot evenness vs. campy for all samples
p = ggplot(copy_sample_data[copy_sample_data$infected=="Infected",], aes(y=log(Observed),x=campy, color=days_post_infection)) + facet_wrap(~diet) +geom_point(size=4,alpha=0.5)
p3 = p + theme(strip.text.x=element_text(size=size_striptext,face="bold"),
              axis.title.x = element_text(size=size_title), 
              axis.text.x = element_text(size=size_axistext, angle=45,hjust=1),
              axis.text.y = element_text(size=size_axistext),
              axis.title.y = element_text(size=size_title),
              legend.position="none",
              plot.margin = unit(c(1,1,1,1), "cm")) + labs(x=xlabel_expression)
ggsave(file="~/Documents/projects/campy_murine_diets/results/campy_v_Observed.svg", plot=p3, width=8, height=4)

# plot evenness vs. campy for all samples
p = ggplot(copy_sample_data[copy_sample_data$infected=="Infected",], aes(y=log(InvSimpson),x=campy, color=days_post_infection)) + facet_wrap(~diet) +geom_point(size=4,alpha=0.5)
p4 = p + theme(strip.text.x=element_text(size=size_striptext,face="bold"),
              axis.title.x = element_text(size=size_title), 
              axis.text.x = element_text(size=size_axistext, angle=45,hjust=1),
              axis.text.y = element_text(size=size_axistext),
              axis.title.y = element_text(size=size_title),
              legend.position="none",
              plot.margin = unit(c(1,1,1,1), "cm")) + labs(x=xlabel_expression) +
              scale_fill_discrete(guide=FALSE)
ggsave(file="~/Documents/projects/campy_murine_diets/results/campy_v_InvSimpson.svg", plot=p4, width=8, height=4)

###------------------------------------------------------------------------

### NMDS of beta diversity using rarefied data
ps_rare
for_factor_conversion = sample_data(ps_rare)
for_factor_conversion$days_post_infection = as.factor(for_factor_conversion$days_post_infection)
sample_data(ps_rare) = for_factor_conversion
ps_rare.ord = ordinate(ps_rare,  method="NMDS", distance = "bray")

p5 = plot_ordination(
  physeq = ps_rare,
  ordination = ps_rare.ord,
  color = "days_post_infection",
  shape = "infected"
) + facet_wrap(~diet) + geom_point(size=6,alpha = 0.5) + 
  theme(strip.text.x=element_text(size=size_striptext,face="bold"),
        axis.title.x = element_text(size=size_title), 
        axis.text.x = element_text(size=size_axistext, angle=45,hjust=1),
        axis.text.y = element_text(size=size_axistext),
        axis.title.y = element_text(size=size_title),
        legend.position="none",
        plot.margin = unit(c(1,1,1,1), "cm"))
###
ggsave(file="~/Documents/projects/campy_murine_diets/results/nmds_all.svg", plot=p5, width=12, height=4)


### Save multipanel plot for Figure 1
pgrid = grid.arrange(p1,p2,p5,nrow=2,
                     layout_matrix=rbind(c(1,2),c(3,3)),
                     heights=c(1.2,1.5))
ggsave(file="~/Documents/projects/campy_murine_diets/results/fig1_multipanel.svg",plot=pgrid,width=16, height=11)

### DIFFERENTIAL ABUNDANCE-------------------------------------------------

# We will use filtered data with low-abundance, rare ASVs removed for differential abundance
rare_min = 5000
# remove samples with fewer counts than the minimum
ps_rare_diff = prune_samples(sample_sums(ps_trimmed) >=rare_min,ps_trimmed)
ps_rare_diff = rarefy_even_depth(ps_rare_diff,sample.size=5000,rngseed=SEED_INT) # OTUs no longer in the dataset will be removed


library("DESeq2")
alpha = 0.05



### Compare each diet at d0 first, since microbiota might be implicated in susceptibility
ps_d0 = prune_samples(sample_data(ps_rare_diff)$days_post_infection %in% c(0),ps_rare_diff)
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
ps_d0 = prune_samples(sample_data(ps_rare_diff)$days_post_infection %in% c(0),ps_rare_diff)
ps_d0_taxatrim = filter_taxa(ps_d0, function(x) sum(x > 1) > 3, TRUE)
ps_d0_taxatrim = filter_taxa(ps_d0_taxatrim, function(x) max(x) > 10, TRUE)
ps_d1 = prune_samples(sample_data(ps_rare_diff)$days_post_infection %in% c(1),ps_rare_diff)
ps_d1_taxatrim = filter_taxa(ps_d1, function(x) sum(x > 1) > 3, TRUE)
ps_d1_taxatrim = filter_taxa(ps_d1_taxatrim, function(x) max(x) > 10, TRUE)
ps_d5 = prune_samples(sample_data(ps_rare_diff)$days_post_infection %in% c(5),ps_rare_diff)
ps_d5_taxatrim = filter_taxa(ps_d5, function(x) sum(x > 1) > 3, TRUE)
ps_d5_taxatrim = filter_taxa(ps_d5_taxatrim, function(x) max(x) > 10, TRUE)
ps_d9 = prune_samples(sample_data(ps_rare_diff)$days_post_infection %in% c(9),ps_rare_diff)
ps_d9_taxatrim = filter_taxa(ps_d9, function(x) sum(x > 1) > 3, TRUE)
ps_d9_taxatrim = filter_taxa(ps_d9_taxatrim, function(x) max(x) > 10, TRUE)

#infection_comparison = function(x,ps_obj,path,alpha) {
#  comp = prune_samples(sample_data(ps_obj)$diet %in% c(x),ps_obj)
#  comp = filter_taxa(comp, function(x) sum(x > 1) > 3, TRUE)
#  
#  comp_deseq = phyloseq_to_deseq2(comp, ~ infected)
#  print(comp_deseq)
#  comp_deseq = DESeq(comp_deseq, test="Wald", fitType="parametric")
#  comp_deseq = results(comp_deseq)
#  print(sum(is.na(comp_deseq)))
#  comp_deseq = merge(comp_deseq,tax_table(comp)[rownames(comp_deseq),],by="row.names")
#  print(sum(is.na(comp_deseq)))
#  # replace NA p values with 1
#  comp_deseq$padj[is.na(comp_deseq$padj)] = 1.0
#  comp_deseq = comp_deseq[comp_deseq$padj < alpha,]
  # save
 # write.table(comp_deseq, file = path, sep = "\t", quote = FALSE,row.names = FALSE)
#  return(comp_deseq)
#}
# run the comparison for all time points and save the results.
#d1_dN = infection_comparison(x="dN",ps_obj=ps_d1_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d1_dN_infection_diff.txt",alpha=0.05)
#d5_dN = infection_comparison(x="dN",ps_obj=ps_d5_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d5_dN_infection_diff.txt",alpha=0.05)
#d9_dN = infection_comparison(x="dN",ps_obj=ps_d9_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d9_dN_infection_diff.txt",alpha=0.05)
#d0_dPD = infection_comparison(x="dPD",ps_obj=ps_d0_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d0_dPD_infection_diff.txt",alpha=0.05)
#d1_dPD = infection_comparison(x="dPD",ps_obj=ps_d1_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d1_dPD_infection_diff.txt",alpha=0.05)
#d5_dPD = infection_comparison(x="dPD",ps_obj=ps_d5_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d5_dPD_infection_diff.txt",alpha=0.05)
#d9_dPD = infection_comparison(x="dPD",ps_obj=ps_d9_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d9_dPD_infection_diff.txt",alpha=0.05)
#d0_dZD = infection_comparison(x="dZD",ps_obj=ps_d0_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d0_dZD_infection_diff.txt",alpha=0.05)
#d1_dZD = infection_comparison(x="dZD",ps_obj=ps_d1_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d1_dZD_infection_diff.txt",alpha=0.05)
#d5_dZD = infection_comparison(x="dZD",ps_obj=ps_d5_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d5_dZD_infection_diff.txt",alpha=0.05)
#d9_dZD = infection_comparison(x="dZD",ps_obj=ps_d9_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d9_dZD_infection_diff.txt",alpha=0.05)
#d0_HC = infection_comparison(x="HC",ps_obj=ps_d0_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d0_HC_infection_diff.txt",alpha=0.05)
#d1_HC = infection_comparison(x="HC",ps_obj=ps_d1_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d1_HC_infection_diff.txt",alpha=0.05)
#d5_HC = infection_comparison(x="HC",ps_obj=ps_d5_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d5_HC_infection_diff.txt",alpha=0.05)
#d9_HC = infection_comparison(x="HC",ps_obj=ps_d9_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d9_HC_infection_diff.txt",alpha=0.05)

# try with a single DESeq2 run
ps_rare_diff_all = ps_rare_diff
sample_data(ps_rare_diff_all)$days_post_infection = as.factor(sample_data(ps_rare_diff_all)$days_post_infection)
sample_data(ps_rare_diff_all)$deseq_group = factor(paste0(sample_data(ps_rare_diff_all)$infected, sample_data(ps_rare_diff_all)$diet,sample_data(ps_rare_diff_all)$days_post_infection))
comp_deseq = phyloseq_to_deseq2(ps_rare_diff_all, ~deseq_group)
comp_deseq = DESeq(comp_deseq, test="Wald", fitType="parametric")
resultsNames(comp_deseq)
#define a function to extract the deseq results
extract_deseq_results = function (comp_deseq,phyloseq_obj,comparisons,alpha) {
  comparison_results = results(comp_deseq,comparisons)
  comparison_frame = merge(comparison_results,tax_table(phyloseq_obj)[rownames(comparison_results),],by="row.names")
  # replace NA p values with 1
  comparison_frame$padj[is.na(comparison_frame$padj)] = 1.0
  comparison_frame = comparison_frame[comparison_frame$padj < alpha,]
}
dZD_0dpi_frame = extract_deseq_results(comp_deseq,ps_rare_diff_all,c("deseq_group","FALSEdZD0","TRUEdZD0"),0.05)
dZD_1dpi_frame = extract_deseq_results(comp_deseq,ps_rare_diff_all,c("deseq_group","FALSEdZD1","TRUEdZD1"),0.05)
dZD_5dpi_frame = extract_deseq_results(comp_deseq,ps_rare_diff_all,c("deseq_group","FALSEdZD5","TRUEdZD5"),0.05)
dZD_9dpi_frame = extract_deseq_results(comp_deseq,ps_rare_diff_all,c("deseq_group","FALSEdZD9","TRUEdZD9"),0.05)

dPD_0dpi_frame = extract_deseq_results(comp_deseq,ps_rare_diff_all,c("deseq_group","FALSEdPD0","TRUEdPD0"),0.05)
dPD_1dpi_frame = extract_deseq_results(comp_deseq,ps_rare_diff_all,c("deseq_group","FALSEdPD1","TRUEdPD1"),0.05)
dPD_5dpi_frame = extract_deseq_results(comp_deseq,ps_rare_diff_all,c("deseq_group","FALSEdPD5","TRUEdPD5"),0.05)
dPD_9dpi_frame = extract_deseq_results(comp_deseq,ps_rare_diff_all,c("deseq_group","FALSEdPD9","TRUEdPD9"),0.05)

HC_0dpi_frame = extract_deseq_results(comp_deseq,ps_rare_diff_all,c("deseq_group","FALSEHC0","TRUEHC0"),0.05)
HC_1dpi_frame = extract_deseq_results(comp_deseq,ps_rare_diff_all,c("deseq_group","FALSEHC1","TRUEHC1"),0.05)
HC_5dpi_frame = extract_deseq_results(comp_deseq,ps_rare_diff_all,c("deseq_group","FALSEHC5","TRUEHC5"),0.05)
HC_9dpi_frame = extract_deseq_results(comp_deseq,ps_rare_diff_all,c("deseq_group","FALSEHC9","TRUEHC9"),0.05)

all_comps_running = merge(dZD_0dpi_frame,dZD_1dpi_frame,c("Kingdom","Phylum","Class","Order","Family","Genus","Row.names"),all=T,suffixes=c('dzd0','dzd1'))
all_taxa = union(dZD_0dpi_frame$Row.names,dZD_1dpi_frame$Row.names)
all_taxa = union(all_taxa,dZD_5dpi_frame$Row.names)
all_taxa = union(all_taxa,dZD_9dpi_frame$Row.names)
all_taxa = union(all_taxa,dPD_0dpi_frame$Row.names)
all_taxa = union(all_taxa,dPD_1dpi_frame$Row.names)
all_taxa = union(all_taxa,dPD_5dpi_frame$Row.names)
all_taxa = union(all_taxa,dPD_9dpi_frame$Row.names)
all_taxa = union(all_taxa,HC_0dpi_frame$Row.names)
all_taxa = union(all_taxa,HC_1dpi_frame$Row.names)
all_taxa = union(all_taxa,HC_5dpi_frame$Row.names)
all_taxa = union(all_taxa,HC_9dpi_frame$Row.names)

alldf = as.data.frame(all_taxa)
alldf$Row.names = alldf$all_taxa
master_results = alldf
ress = merge(dZD_0dpi_frame,alldf,"Row.names",all=T)
master_results$dzd0_fc = ress$log2FoldChange
ress = merge(dZD_1dpi_frame,alldf,"Row.names",all=T)
master_results$dzd1_fc = ress$log2FoldChange
ress = merge(dZD_5dpi_frame,alldf,"Row.names",all=T)
master_results$dzd5_fc = ress$log2FoldChange
ress = merge(dZD_9dpi_frame,alldf,"Row.names",all=T)
master_results$dzd9_fc = ress$log2FoldChange

ress = merge(dPD_0dpi_frame,alldf,"Row.names",all=T)
master_results$dpd0_fc = ress$log2FoldChange
ress = merge(dPD_1dpi_frame,alldf,"Row.names",all=T)
master_results$dpd1_fc = ress$log2FoldChange
ress = merge(dPD_5dpi_frame,alldf,"Row.names",all=T)
master_results$dpd5_fc = ress$log2FoldChange
ress = merge(dPD_9dpi_frame,alldf,"Row.names",all=T)
master_results$dpd9_fc = ress$log2FoldChange

ress = merge(HC_0dpi_frame,alldf,"Row.names",all=T)
master_results$hc0_fc = ress$log2FoldChange
ress = merge(HC_1dpi_frame,alldf,"Row.names",all=T)
master_results$hc1_fc = ress$log2FoldChange
ress = merge(HC_5dpi_frame,alldf,"Row.names",all=T)
master_results$hc5_fc = ress$log2FoldChange
ress = merge(HC_9dpi_frame,alldf,"Row.names",all=T)
master_results$hc9_fc = ress$log2FoldChange
# replace NA with 0
master_results[is.na(master_results)] = 0

row.names(master_results) = paste0(tax_table(ps_rare_diff_all)[master_results$Row.names,'Family'],
                                   ' ',tax_table(ps_rare_diff_all)[master_results$Row.names,'Genus'],
                                   ' ',nrow(master_results)+1-seq_along(taxa_names(ps_rare_diff_all)[1:nrow(master_results)]))
master_results$all_taxa = NULL
master_results$Row.names = NULL
paletteLength = 100
myBreaks <- c(seq(min(master_results), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(master_results)/paletteLength, max(master_results), length.out=floor(paletteLength/2)))
pheatmap(master_results,cluster_cols=FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(paletteLength),
         breaks = myBreaks,
         file="~/Documents/projects/campy_murine_diets/results/diffabund_heatmap.png",
         width = 8,
         height = 6)

