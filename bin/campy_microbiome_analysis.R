library(dada2)
library(ShortRead)
library(ggplot2)
library(phyloseq)
library(reshape2)


ps <- readRDS("~/Documents/projects/campy_murine_diets/data/campy_phyloseq_obj_less_strict.rds")
neg_control = prune_samples(rownames(sample_data(ps)) %in% c("Plate3-D03"),ps)
path = "~/Documents/projects/campy_murine_diets/data/"
metadata = read.table(paste0(path,"Campy3microbiome1-23-2017.txt"),header=TRUE,
           sep = "\t")
colnames(metadata)[colnames(metadata) == 'Plate'] = "plate_sample"
weights = read.table(paste0(path,"campystudy_weights.txt"),header=TRUE,
                      sep = "\t",check.names=FALSE)
#shedding = read.table(paste0(path,"campystudy_shedding.txt"),header=TRUE,
#                     sep = "\t",check.names = FALSE)
#shedding_long = melt(shedding,id.vars =c('Day'))
#colnames(shedding_long)[colnames(shedding_long) == "variable"] = "Sample"
#colnames(shedding_long)[colnames(shedding_long) == "value"] = "shedding"
#colnames(shedding_long)[colnames(shedding_long) == "Day"] = "days_post_infection"

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
# subset the amples in the phyloseq object, then reassign the metadata as sample_data
ps = prune_samples(rownames(sample_data(ps)) %in% rownames(new_sampledata), ps)
#sample_data(ps) = new_sampledata
sample_data(ps) = new_sampledata

# prune samples with fewer than 10 counts, and those with fewer than 5 non-zero OTUS
old_samples = rownames(sample_data(ps))
ps = prune_samples(sample_sums(ps)>=10, ps)
# get samples with fewer than 5 non-zero OTUs
nonzero_count = apply(otu_table(ps), 1, function(c)sum(c!=0))

ps = prune_samples(nonzero_count > 5, ps)

pruned_samples = old_samples[!old_samples %in% rownames(sample_data(ps))]

old_taxa = rownames(tax_table(ps))
ps_count_taxatrim = filter_taxa(ps, function(x) max(x) > 20, TRUE)
ps_count_taxatrim <- filter_taxa(ps_count_taxatrim, function(x) sum(x > 1) > 5, TRUE)


ps_relative_taxatrim  = transform_sample_counts(ps_count_taxatrim, function(x) x / sum(x) )

### figure 1: alpha and beta diversity.
# Figure 1a: alpha diversity
ps_d015 = prune_samples(sample_data(ps)$days_post_infection %in% c(0,1,5),ps)
# prune the taxa that don't appear in any of these samples (e.g. only found in day 9&11)
ps_d015 <- filter_taxa(ps_d015, function(x) sum(x > 1) > 1, TRUE)

sample_data(ps_d015)$days_post_infection = as.factor(sample_data(ps_d015)$days_post_infection)
# get actual values and test differences
richness = estimate_richness(ps_d015,measures="Shannon")
asdf = as(sample_data(ps_d015),"data.frame")
asdf = merge(asdf,richness,by="row.names")
asdf$days_post_infection = as.factor(asdf$days_post_infection)
asdf$infected = as.factor(asdf$infected)
diet_anova = aov(Shannon~diet*days_post_infection*infected,data=asdf)
diet_tukey_posthoc = TukeyHSD(diet_anova)

theme_set(theme_bw())
p = plot_richness(ps_d015,measures = "Shannon", x='days_post_infection')
p = p + geom_boxplot(alpha=0.5) + facet_grid(infected~diet)
#p$layers = p$layers[-1]
p
# save as svg
ggsave(file="~/Documents/projects/campy_murine_diets/results/Shannon_div_all.svg", plot=p, width=10, height=8)
# save the anova results
saveRDS(diet_tukey_posthoc, "~/Documents/projects/campy_murine_diets/results/Shannon_anova_posthoc.rds")


# Figure 1b: beta diversity of mice pre and post-infection
# 
ps_d015.ord <- ordinate(ps_d015, "NMDS", "bray")
theme_set(theme_bw())
p1 = plot_ordination(ps_d015, ps_d015.ord, color="days_post_infection",type="samples", shape="infected") + facet_wrap(~diet) + theme(strip.text = element_text(size=20)) 
p1 = p1 + geom_point(size=3,alpha=1)
#new_data = newdata=data.frame(ps_d015.ord$points[,1],ps_d015.ord$points[,2],sample_data(ps_d015)$infected)
#colnames(new_data) = c('NMDS1','NMDS2','infected')
#p1 = p1 + geom_point(data=new_data,mapping=aes(shape='infected'))
p1
ggsave(file="~/Documents/projects/campy_murine_diets/results/beta_div_bray_all.svg", plot=p1, width=10, height=8)

# Differential abundance: infected vs. uninfected @ each time point
library("DESeq2")
alpha = 0.05
# trim low abundance taxa and samples. Must be in 5 samples, at least max of 10 counts in one sample.
ps_d015_taxatrim = filter_taxa(ps_d015, function(x) sum(x > 1) > 5, TRUE)
ps_d015_taxatrim = filter_taxa(ps_d015_taxatrim, function(x) max(x) > 10, TRUE)

# Compare each diet at d0
ps_d0 = prune_samples(sample_data(ps_d015)$days_post_infection %in% c(0),ps_d015)
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
d0_dNvdPD = pair_diet_comparison(x="dN",y="dPD",path="~/Documents/projects/campy_murine_diets/results/d0_dNvdPD_diff.txt",alpha=0.05)
d0_dNvdZD = pair_diet_comparison(x="dN",y="dZD",path="~/Documents/projects/campy_murine_diets/results/d0_dNvdZD_diff.txt",alpha=0.05)
d0_dNvHC = pair_diet_comparison(x="dN",y="HC",path="~/Documents/projects/campy_murine_diets/results/d0_dNvHC_diff.txt",alpha=0.05)
d0_dZDvdPD = pair_diet_comparison(x="dZD",y="dPD",path="~/Documents/projects/campy_murine_diets/results/d0_dZDvdPD_diff.txt",alpha=0.05)
d0_dZDvHC = pair_diet_comparison(x="dZD",y="HC",path="~/Documents/projects/campy_murine_diets/results/d0_dZDvHC_diff.txt",alpha=0.05)
d0_dPDvHC = pair_diet_comparison(x="dPD",y="HC",path="~/Documents/projects/campy_murine_diets/results/d0_dPDvHC_diff.txt",alpha=0.05)

# difference between dZD and the rest?
higher_in_dZD = d0_dZDvdPD[d0_dZDvdPD$log2FoldChange > 0,]
higher_in_dZD = d0_dZDvHC[d0_dZDvHC$log2FoldChange > 0,]
higher_in_dZD = d0_dNvdZD[d0_dNvdZD$log2FoldChange < 0,]

lower_in_dZD = d0_dZDvdPD[d0_dZDvdPD$log2FoldChange < 0,][c("Row.names","Order","Family","Genus")]
lower_in_dZD = merge(lower_in_dZD,d0_dZDvHC[d0_dZDvHC$log2FoldChange < 0,][c("Row.names","Order","Family","Genus")],by="Row.names")
lower_in_dZD = merge(lower_in_dZD,d0_dNvdZD[d0_dNvdZD$log2FoldChange > 0,][c("Row.names","Order","Family","Genus")],by="Row.names")

# Effect of infection at each time point, for each diet
ps_d1 = prune_samples(sample_data(ps_d015)$days_post_infection %in% c(1),ps_d015)
ps_d1_taxatrim = filter_taxa(ps_d1, function(x) sum(x > 1) > 3, TRUE)
ps_d1_taxatrim = filter_taxa(ps_d1_taxatrim, function(x) max(x) > 10, TRUE)
ps_d5 = prune_samples(sample_data(ps_d015)$days_post_infection %in% c(5),ps_d015)
ps_d5_taxatrim = filter_taxa(ps_d5, function(x) sum(x > 1) > 3, TRUE)
ps_d5_taxatrim = filter_taxa(ps_d5_taxatrim, function(x) max(x) > 10, TRUE)

infection_comparison = function(x,ps_obj,path,alpha) {
  comp = prune_samples(sample_data(ps_obj)$diet %in% c(x),ps_obj)
  comp = filter_taxa(comp, function(x) sum(x > 1) > 3, TRUE)
  comp_deseq = phyloseq_to_deseq2(comp, ~ infected)
  comp_deseq = DESeq(comp_deseq, test="Wald", fitType="parametric")
  comp_deseq = results(comp_deseq)
  #print(dim(comp_deseq))
  comp_deseq = merge(comp_deseq,tax_table(comp)[rownames(comp_deseq),],by="row.names")
  #print(dim(comp_deseq))
  comp_deseq = comp_deseq[comp_deseq$padj < alpha,]
  # save
  write.table(comp_deseq, file = path, sep = "\t", quote = FALSE,row.names = FALSE)
  return(comp_deseq)
}
d1_dN = infection_comparison(x="dN",ps_obj=ps_d1_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d1_dN_infection_diff.txt",alpha=0.05)
d5_dN = infection_comparison(x="dN",ps_obj=ps_d5_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d5_dN_infection_diff.txt",alpha=0.05)
d1_dPD = infection_comparison(x="dPD",ps_obj=ps_d1_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d1_dPD_infection_diff.txt",alpha=0.05)
d5_dPD = infection_comparison(x="dPD",ps_obj=ps_d5_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d5_dPD_infection_diff.txt",alpha=0.05)
d1_dZD = infection_comparison(x="dZD",ps_obj=ps_d1_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d1_dZD_infection_diff.txt",alpha=0.05)
d5_dZD = infection_comparison(x="dZD",ps_obj=ps_d5_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d5_dZD_infection_diff.txt",alpha=0.05)
d1_HC = infection_comparison(x="HC",ps_obj=ps_d1_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d1_HC_infection_diff.txt",alpha=0.05)
d5_HC = infection_comparison(x="HC",ps_obj=ps_d5_taxatrim,path="~/Documents/projects/campy_murine_diets/results/d5_HC_infection_diff.txt",alpha=0.05)


## Save OTU table, sample data, and taxonomy table for natasa
ps_as_matrix = as(otu_table(ps_relative_taxatrim),"matrix")
ps_as_df = as.data.frame(ps_as_matrix)
taxa_as_matrix = as(tax_table(ps_relative_taxatrim),"matrix")
taxa_as_df = as.data.frame(taxa_as_matrix)
sampledata_as_matrix = as(sample_data(ps_relative_taxatrim),"matrix")
sampledata_as_df = as.data.frame(sampledata_as_matrix)
# rename the columns of relative abundance df with family_genus_num
colnames(ps_as_df) = paste0(tax_table(ps_relative_taxatrim)[,'Family'],tax_table(ps_relative_taxatrim)[,'Genus'],seq_along(taxa_names(ps_relative_taxatrim)))
taxa_as_df$name_in_abundance_table = colnames(ps_as_df)
# save dataframes
write.table(ps_as_df,file='~/Documents/projects/campy_murine_diets/results/abundance_table.tsv',sep='\t',col.names=NA)
write.table(taxa_as_df,file='~/Documents/projects/campy_murine_diets/results/taxa_table.tsv',sep='\t',col.names=NA)
write.table(sampledata_as_df,file='~/Documents/projects/campy_murine_diets/results/sampledata_table.tsv',sep='\t',col.names=NA)

#ps_d0_taxatrim_deseq = phyloseq_to_deseq2(ps_d0_taxatrim, ~ diet)
#ps_d0_taxatrim_deseq = DESeq(ps_d0_taxatrim_deseq, test="Wald", fitType="parametric")
#ps_d0_taxatrim_deseq = results(ps_d0_taxatrim_deseq)


#ps_d015_taxatrim_deseq = phyloseq_to_deseq2(ps_d015_taxatrim, ~ infected + days_post_infection)
#ps_d015_taxatrim_deseq = DESeq(ps_d015_taxatrim_deseq, test="Wald", fitType="parametric")
#ps_d015_taxatrim_deseq = results(ps_d015_taxatrim_deseq)
#dZD_d1_inf_sigtab = dZD_d1_inf_deseq[which(dZD_d1_inf_deseq$padj < alpha), ]
#dZD_d1_inf_sigtab = cbind(as(dZD_d1_inf_sigtab, "data.frame"), as(tax_table(ps_count_taxatrim_dZD_d1)[rownames(dZD_d1_inf_sigtab), ], "matrix"))




#plot just dZD
#ps_relative_taxatrim_dZD = prune_samples(sample_data(ps_relative_taxatrim)$diet == 'dZD', ps_relative_taxatrim)
#ps_relative_taxatrim_dPD = prune_samples(sample_data(ps_relative_taxatrim)$diet == 'dPD', ps_relative_taxatrim)
#ps_relative_taxatrim_dN = prune_samples(sample_data(ps_relative_taxatrim)$diet == 'dN', ps_relative_taxatrim)
#ps_relative_taxatrim_HC = prune_samples(sample_data(ps_relative_taxatrim)$diet == 'HC', ps_relative_taxatrim)


#plot_bar(ps_relative_taxatrim_dZD, x="days_post_infection", fill="Phylum") + facet_grid(infected~diet, scales="free_x")

#campy_seq = rownames(tax_table(ps)[which(tax_table(ps)[,'Genus'] == 'Campylobacter')])
#otu_table(ps)[,colnames(otu_table(ps))==campy_seq]

#alt_campy_seqs = rownames(tax_table(ps)[which(tax_table(ps)[,"Phylum"] == "Proteobacteria")])
#otu_table(ps)[,colnames(otu_table(ps))==alt_campy_seqs[2]]

# check negative control
# filter by non-zero taxa
neg_control = filter_taxa(neg_control, function(x) x > 0, TRUE)
neg_control_table = tax_table(neg_control)
View(neg_control_table)



# rename taxa for heatmap plotting
taxa_names(ps_count_taxatrim_d0) = paste0(tax_table(ps_count_taxatrim_d0)[,'Family'],tax_table(ps_count_taxatrim_d0)[,'Genus'],seq_along(taxa_names(ps_count_taxatrim_d0)))
labels = sample_data(ps_count_taxatrim_d0)
plot_heatmap(ps_count_taxatrim_d0,sample.label = "diet")

# generate predictive model for diet to learn dZD predictive seqs.
ps_relative_taxatrim_shedding = prune_samples(!is.na(sample_data(ps_relative_taxatrim)$campy),ps_relative_taxatrim)
# remove the samples from uninfected groups, and stop at day 5
only_infected = prune_samples(sample_data(ps_relative_taxatrim_shedding)$infected == TRUE,ps_relative_taxatrim_shedding)
only_infected = prune_samples(sample_data(only_infected)$days_post_infection <= 5,only_infected)
only_infected = prune_samples(sample_data(only_infected)$days_post_infection == 1,only_infected)
only_infected = prune_samples(sample_data(only_infected)$diet %in% c('dPD','dN','dZD'),only_infected)

for_rf = prune_samples(sample_data(ps_relative_taxatrim)$days_post_infection == 5,ps_relative_taxatrim_shedding)
x = as(otu_table(for_rf),"matrix")
#x = cbind(x,sample_data(only_infected)$days_post_infection)
y = sample_data(for_rf)$campy
# non-detects were entered as 10, reassign to 0.
y[y < 50] = FALSE
y[y > 50] = TRUE
y = as.factor(y)
# generate predictive model for shedding and weights
##weight_y = sample_data(only_infected)$campy


library("randomForest")
rf = randomForest(x=x,y=y)
#weight_rf = randomForest(x=x,y=weight_y)
row.names(rf$importance) = paste0(tax_table(for_rf)[,'Family'][row.names(rf$importance)],tax_table(for_rf)[,'Genus'][row.names(rf$importance)],seq_along(taxa_names(for_rf)))
#row.names(weight_rf$importance) = paste0(tax_table(only_infected)[,'Family'][row.names(weight_rf$importance)],tax_table(only_infected)[,'Genus'][row.names(weight_rf$importance)],seq_along(taxa_names(only_infected)))

# plot heatmap using only highly important OTUs
important = rf$importance
top10 = order(rf$importance,decreasing = TRUE)[1:10]
top10_names = rownames(rf$importance)[order(rf$importance, decreasing=TRUE)][1:10]

ps_filt_by_imp = prune_taxa(taxa = top10_names,for_rf)
taxa_names(ps_filt_by_imp) = paste0(tax_table(for_rf)[,'Family'][top10_names],tax_table(only_infected)[,'Genus'][top10_names],seq_along(taxa_names(ps_filt_by_imp)))
plot_heatmap(ps_filt_by_imp,sample.label = "diet")
partialPlot(rf,pred.data = x ,x.var = top10[1])
# for each of the RF-important SVs, make a scatterplot by diet, time, and infected/uninfected



# plot lactobacillus aboundance in all samples
# get lactobacillus abundance and add to sample data dataframe
meta = sample_data(for_rf)
# get important seq values and add to meta
important_seqs = as(otu_table(for_rf),"matrix")[,top10]
# rename columns with the taxonomy
colnames(important_seqs) = paste0(tax_table(for_rf)[,'Family'][top10],tax_table(for_rf)[,'Genus'][top10],seq_along(colnames(important_seqs)))
#lacto = as(otu_table(ps_relative_taxatrim)[,lacto_seq],"matrix") # samples are rows
#meta$lacto_abundance = lacto
meta[colnames(important_seqs)] = important_seqs
meta = as.data.frame(meta)
meta$days_post_infection = as.factor(meta$days_post_infection)

# plot scatter by time, diet, and infection
#enterococcus 8, lactobacillus 15
seq = colnames(important_seqs)[1]
p <- ggplot(meta, aes_string(x="days_post_infection",y=seq)) + geom_point() + facet_grid(~diet, scales="free_x")
p <- ggplot(meta, aes_string(x="days_post_infection",y='campy',colour='infected')) + geom_point() + facet_grid(~diet, scales="free_x")
p
meta_copy = meta
meta_copy$infected = as.factor(meta_copy$infected)
plot(meta_copy$EnterobacteriaceaeNA1, meta_copy$EnterobacteriaceaeNA2,col=meta_copy$infected)

p <- ggplot(meta_copy[meta_copy$infected==TRUE,], aes(x=campy, y=ErysipelotrichaceaeTuricibacter4,colour=diet)) + geom_point()
p <- ggplot(meta_copy[meta_copy$infected==TRUE,], aes(x=campy, y=Bacteroidales_S24_7_groupNA6,colour=diet)) + geom_point()
p <- ggplot(meta_copy[meta_copy$infected==TRUE,], aes(x=campy, y=LachnospiraceaeMarvinbryantia7,colour=diet)) + geom_point()
p <- ggplot(meta_copy[meta_copy$infected==TRUE,], aes(x=campy, y=LachnospiraceaeNA3,colour=diet)) + geom_point()
p <- ggplot(meta_copy[meta_copy$infected==TRUE,], aes(x=campy, y=LactobacillaceaeLactobacillus2,colour=diet)) + geom_point()
p <- ggplot(meta_copy[meta_copy$infected==TRUE,], aes(x=campy, y=StreptococcaceaeStreptococcus1,colour=diet)) + geom_point()
p <- ggplot(meta_copy, aes(x=campy, y=VerrucomicrobiaceaeAkkermansia9,colour=diet)) + geom_point()# + facet_grid(~diet, scales="free_x")
p <- ggplot(meta_copy, aes(x=campy, y=Clostridiales_vadinBB60_groupNA1,colour=diet)) + geom_point()# + facet_grid(~diet, scales="free_x")

p

