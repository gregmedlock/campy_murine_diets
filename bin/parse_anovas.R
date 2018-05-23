
diet_tukey_posthoc_richness = readRDS("~/Documents/projects/campy_murine_diets/results/Observed_anova_posthoc.rds")
diet_tukey_posthoc_evenness = readRDS("~/Documents/projects/campy_murine_diets/results/InvSimpson_anova_posthoc.rds")

threeway_richness = data.frame(diet_tukey_posthoc_richness$`diet:days_post_infection:infected`)
threeway_evenness = data.frame(diet_tukey_posthoc_evenness$`diet:days_post_infection:infected`)

threeway_richness_sig = threeway_richness[threeway_richness$p.adj < 0.05,]
threeway_evenness_sig = threeway_evenness[threeway_evenness$p.adj < 0.05,]


time_diffs = function(alpha_data,diets_excluded) {
  alpha_data[!grepl(paste(diets_excluded,collapse="|"),rownames(alpha_data)),]
}

dZD_time_evenness = time_diffs(threeway_evenness_sig,c("dPD","HC"))
dPD_time_evenness = time_diffs(threeway_evenness_sig,c("dZD","HC"))
HC_time_evenness = time_diffs(threeway_evenness_sig,c("dPD","dZD"))

dZD_time_richness = time_diffs(threeway_richness_sig,c("dPD","HC"))
dPD_time_richness = time_diffs(threeway_richness_sig,c("dZD","HC"))
HC_time_richness = time_diffs(threeway_richness_sig,c("dPD","dZD"))
