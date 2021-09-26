

samples_lungmetastasis = readtable('/Users/s202425/Documents/GitHub/metastasis/results/sampling/sampling_lung_metastasis.csv');
samples_brainmetastasis = readtable('/Users/s202425/Documents/GitHub/metastasis/results/sampling/sampling_brain_metastasis.csv');
samples_breastcancer = readtable('/Users/s202425/Documents/GitHub/metastasis/results/sampling/sampling_breast_cancer.csv');


% LUNGEMETASTASIS VS BRAINMETASTASIS
% find common reactions
[commonreaction, pos1] = intersect(samples_lungmetastasis.Properties.VariableNames, samples_brainmetastasis.Properties.VariableNames)

file = fopen('/Users/s202425/Documents/GitHub/metastasis/results/sampling_statisticalcomparison/lungmetastasis_vs_brainmetastasis.csv', 'a');
fprintf(file, 'reaction;jaccardindex;KolmogorovSmirnow;KolmogorovSmirenow_pvalue;T-test_pvalue;Anova_pvalue;Wilcoxon_Mann_Whitney_pvalue;Mean;Std;Mode;Median;Skew;minFlux_lungmetastasis;minFlux_brainmetastasis;maxFlux_lungmetastasis;maxFlux_brainmetastasis\n');

for pos = 1:length(commonreaction)
   reaction = commonreaction{pos};
   % get sample array from specific array
   samp_lungmetastasis = samples_lungmetastasis{:, {reaction}};
   samp_brainmetastasis = samples_brainmetastasis{:, {reaction}};
   % Kolmogorov-Smirnov test, rank-sum test, T-test
   [stats, pVals] = compareTwoSamplesStat(samp_lungmetastasis,samp_brainmetastasis, ...
   {'ks', 'rankSum','tTest'});
   sampleStats = calcSampleStats({samp_lungmetastasis, samp_brainmetastasis});
   
   % Wilcoxon Mann Whitney Test
   p_ranksum = ranksum(samp_lungmetastasis, samp_brainmetastasis);
   % jaccard index
   minFlux_lungmet = prctile(samp_lungmetastasis, 25);
   minFlux_brainmetastasis = prctile(samp_brainmetastasis, 25);
   maxFlux_lungmet = prctile(samp_lungmetastasis, 75);
   maxFlux_brainmetastasis = prctile(samp_brainmetastasis, 75);
   minFlux =[minFlux_lungmet(:), minFlux_brainmetastasis(:)];
   maxFlux =[maxFlux_lungmet(:), maxFlux_brainmetastasis(:)];
   
   jaccard_index = fvaJaccardIndex(minFlux, maxFlux);
   % T-test
   [h_ttest, pvalue_ttest]=ttest(samp_lungmetastasis,samp_brainmetastasis);
   % Anova
   matrix_anova = [samp_lungmetastasis(:),samp_brainmetastasis(:)];
   [p_anova, tbl, stats_a] = anova1(matrix_anova);
   filepath = '/Users/s202425/Documents/GitHub/metastasis/results/sampling_statisticalcomparison/boxplot/'  + string(reaction);
   savefig(filepath);
   close all;
   
  
   
   stat_results = [string(jaccard_index),string(mean(stats.ks)),string(mean(pVals.ks)),...
       string(pvalue_ttest), string(p_anova), string(p_ranksum), ...
       strjoin(string(sampleStats.mean), '__'),  strjoin(string(sampleStats.std), '__'),...
       strjoin(string(sampleStats.mode), '__'),strjoin(string(sampleStats.median), '__'),...
       strjoin(string(sampleStats.skew), '__'),strjoin(string(sampleStats.kurt), '__'), ...
       string(h_ttest), string(pvalue_ttest),string(minFlux_lungmet),string(minFlux_brainmetastasis),...
       string(maxFlux_lungmet),string(maxFlux_brainmetastasis)];
   stat_results = fillmissing(stat_results, 'constant', "none");
   
   rowtext = string(reaction) + ';' + strjoin(string(stat_results), ';')+ ';' + '\n';
   fprintf(file, rowtext);
end 
fclose(file);


% LUNGEMETASTASIS VS BREASTCANCER
% find common reactions
[commonreaction, pos1] = intersect(samples_lungmetastasis.Properties.VariableNames, samples_breastcancer.Properties.VariableNames)

file = fopen('/Users/s202425/Documents/GitHub/metastasis/results/sampling_statisticalcomparison/lungmetastasis_vs_breastcancer.csv', 'a');
fprintf(file, 'reaction;jaccardindex;KolmogorovSmirnow;KolmogorovSmirenow_pvalue;T-test_pvalue;Anova_pvalue;Wilcoxon_Mann_Whitney_pvalue;Mean;Std;Mode;Median;Skew;minFlux_lungmetastasis;minFlux_breastcancer;maxFlux_lungmetastasis;maxFlux_breastcancer;perc_lungmet/breastcancer;perc_breastcancer/lungmet\n');

for pos = 1:length(commonreaction)
   reaction = commonreaction{pos};
   % get sample array from specific array
   samp_lungmetastasis = samples_lungmetastasis{:, {reaction}};
   samp_breastcancer = samples_breastcancer{:, {reaction}};
   % Kolmogorov-Smirnov test, rank-sum test, T-test
   [stats, pVals] = compareTwoSamplesStat(samp_lungmetastasis,samp_breastcancer, ...
   {'ks', 'rankSum','tTest'});
   sampleStats = calcSampleStats({samp_lungmetastasis, samp_breastcancer});
   
   % Wilcoxon Mann Whitney Test
   p_ranksum = ranksum(samp_lungmetastasis, samp_breastcancer);
   % jaccard index
   minFlux_lungmet = prctile(samp_lungmetastasis, 25);
   minFlux_breastcancer = prctile(samp_breastcancer, 25);
   maxFlux_lungmet = prctile(samp_lungmetastasis, 75);
   maxFlux_breastcancer = prctile(samp_breastcancer, 75);
   minFlux =[minFlux_lungmet(:), minFlux_breastcancer(:)];
   maxFlux =[maxFlux_lungmet(:), maxFlux_breastcancer(:)];
   
   jaccard_index = fvaJaccardIndex(minFlux, maxFlux);
   % T-test
   [h_ttest, pvalue_ttest]=ttest(samp_lungmetastasis,samp_breastcancer);
   % Anova
   matrix_anova = [samp_lungmetastasis(:),samp_breastcancer(:)];
   [p_anova, tbl, stats_a] = anova1(matrix_anova);
   filepath = '/Users/s202425/Documents/GitHub/metastasis/results/sampling_statisticalcomparison/boxplot_lungmet_breastcancer/'  + string(reaction);
   savefig(filepath);
   close all;

   perc1 = sampleStats.median(1)/sampleStats.median(2);
   perc2 = sampleStats.median(2)/sampleStats.median(1);

   stat_results = [string(jaccard_index),string(mean(stats.ks)),string(mean(pVals.ks)),...
       string(pvalue_ttest), string(p_anova), string(p_ranksum), ...
       strjoin(string(sampleStats.mean), '__'),  strjoin(string(sampleStats.std), '__'),...
       strjoin(string(sampleStats.mode), '__'),strjoin(string(sampleStats.median), '__'),...
       strjoin(string(sampleStats.skew), '__'),strjoin(string(sampleStats.kurt), '__'), ...
       string(h_ttest), string(pvalue_ttest),string(minFlux_lungmet),string(minFlux_breastcancer),...
       string(maxFlux_lungmet),string(maxFlux_breastcancer), perc1, perc2];
   stat_results = fillmissing(stat_results, 'constant', "none");
   
   rowtext = string(reaction) + ';' + strjoin(string(stat_results), ';')+ ';' + '\n';
   fprintf(file, rowtext);
end 
fclose(file);


% BRAINMETASTASIS VS BREASTCANCER
% find common reactions
[commonreaction, pos1] = intersect(samples_brainmetastasis.Properties.VariableNames, samples_breastcancer.Properties.VariableNames)

file = fopen('/Users/s202425/Documents/GitHub/metastasis/results/sampling_statisticalcomparison/brainmetastasis_vs_breastcancer.csv', 'a');
fprintf(file, 'reaction;jaccardindex;KolmogorovSmirnow;KolmogorovSmirenow_pvalue;T-test_pvalue;Anova_pvalue;Wilcoxon_Mann_Whitney_pvalue;Mean;Std;Mode;Median;Skew;minFlux_brainmetastasis;minFlux_breastcancer;maxFlux_brainmetastasis;maxFlux_breastcancer;perc_brainmet/breastcancer;perc_breastcancer/brainmet\n');

for pos = 1:length(commonreaction)
   reaction = commonreaction{pos};
   % get sample array from specific array
   samp_brainmetastasis = samples_brainmetastasis{:, {reaction}};
   samp_breastcancer = samples_breastcancer{:, {reaction}};
   % Kolmogorov-Smirnov test, rank-sum test, T-test
   [stats, pVals] = compareTwoSamplesStat(samp_brainmetastasis,samp_breastcancer, ...
   {'ks', 'rankSum','tTest'});
   sampleStats = calcSampleStats({samp_brainmetastasis, samp_breastcancer});
   
   % Wilcoxon Mann Whitney Test
   p_ranksum = ranksum(samp_brainmetastasis, samp_breastcancer);
   % jaccard index
   minFlux_brainmet = prctile(samp_brainmetastasis, 25);
   minFlux_breastcancer = prctile(samp_breastcancer, 25);
   maxFlux_brainmet = prctile(samp_brainmetastasis, 75);
   maxFlux_breastcancer = prctile(samp_breastcancer, 75);
   minFlux =[minFlux_brainmet(:), minFlux_breastcancer(:)];
   maxFlux =[maxFlux_brainmet(:), maxFlux_breastcancer(:)];
   
   jaccard_index = fvaJaccardIndex(minFlux, maxFlux);
   % T-test
   [h_ttest, pvalue_ttest]=ttest(samp_brainmetastasis,samp_breastcancer);
   % Anova
   matrix_anova = [samp_brainmetastasis(:),samp_breastcancer(:)];
   [p_anova, tbl, stats_a] = anova1(matrix_anova);
   filepath = '/Users/s202425/Documents/GitHub/metastasis/results/sampling_statisticalcomparison/boxplot_brainmet_breastcancer/'  + string(reaction);
   savefig(filepath);
   close all;

   perc1 = sampleStats.median(1)/sampleStats.median(2);
   perc2 = sampleStats.median(2)/sampleStats.median(1);

   stat_results = [string(jaccard_index),string(mean(stats.ks)),string(mean(pVals.ks)),...
       string(pvalue_ttest), string(p_anova), string(p_ranksum), ...
       strjoin(string(sampleStats.mean), '__'),  strjoin(string(sampleStats.std), '__'),...
       strjoin(string(sampleStats.mode), '__'),strjoin(string(sampleStats.median), '__'),...
       strjoin(string(sampleStats.skew), '__'),strjoin(string(sampleStats.kurt), '__'), ...
       string(h_ttest), string(pvalue_ttest),string(minFlux_brainmet),string(minFlux_breastcancer),...
       string(maxFlux_brainmet),string(maxFlux_breastcancer), perc1, perc2];
   stat_results = fillmissing(stat_results, 'constant', "none");
   
   rowtext = string(reaction) + ';' + strjoin(string(stat_results), ';')+ ';' + '\n';
   fprintf(file, rowtext);
end 
fclose(file);




