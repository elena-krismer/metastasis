
[commonreaction, pos1] = intersect(lungtissue.rxns,lung_metastasis.rxns);

samples_lungmetastasis = readtable('/Users/s202425/Documents/GitHub/metastasis/results/sampling/sampling_lung_metastasis.csv');
samples_brainmetastasis = readtable('/Users/s202425/Documents/GitHub/metastasis/results/sampling/sampling_brain_metastasis.csv');



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

x = strjoin(string(sampleStats.mean), ',') + ';' + ...
   strjoin(string(sampleStats.std), ',')

rowtext = reaction + ';' + string(mean(stats.ks)) + ';' + string(mean(stats.ranksum))+';'+ ...
   string(mean(stats.ttest)) + ';' + string(mean(pVals.ks))+ ';'+ string(mean(pVals.ranksum))+ ';'+ ...
   
row = mean(pVals.ttest) + ';' +  strjoin(string(sampleStats.mean), ',') + ';' + ...
   strjoin(string(sampleStats.std), ',') + ...
   ';' + strjoin(string(sampleStats.mode), ',') + ';' + strjoin(string(sampleStats.median), ',') + ';' + ...
   strjoin(string(sampleStats.skew), ',') + ';' + strjoin(string(sampleStats.kurt), ',') + '\n';

% BREASTCANCER VS LUNGTISSUE
%find commonreactions
[commonreaction, pos1] = intersect(samples_lungtissue.Properties.VariableNames, samples_breastcancer.Properties.VariableNames)

file = fopen('/Users/s202425/Documents/GitHub/metastasis/results/sampling/breastcancer_vs_lungtissue.csv', 'a');
fprintf(file, 'reaction;stats.ks;stats.ranksum;stats.ttest;pVals.ks;pVals.ranksum;pVals.ttest;sampleStats.mean;sampleStats.std;sampleStats.mode;sampleStats.median;sampleStats.skew;sampleStats.kurt');

for pos = 1:length(commonreaction)
   reaction = commonreaction{pos};
   % get sample array from specific array 
   samp_breastcancer = samples_breastcancer{:, {reaction}};
   samp_lungtissue = samples_lungtissue{:, {reaction}};
   % Kolmogorov-Smirnov test, rank-sum test, T-test
   [stats, pVals] = compareTwoSamplesStat(samp_breastcancer,samp_lungtissue, {'ks', 'rankSum','tTest'});
   sampleStats = calcSampleStats({samp_breastcancer, samp_lungtissue});
   % write to file
   fprintf(file, reaction, ';', mean(stats.ks), ';', mean(stats.ranksum),';', ...
   mean(stats.ttest), ';', mean(pVals.ks), ';', mean(pVals.ranksum), ';', mean(pVals.ttest), ';',...
   sampleStats.mean, ';', sampleStats.std, ';', sampleStats.mode, ';',...
   sampleStats.median, ';',sampleStats.skew, ';',sampleStats.kurt);
end 
fclose(file);

% BREASTCANCER VS LUNGMETASTASIS
% find common reactions
[commonreaction, pos1] = intersect(samples_lungmet.Properties.VariableNames, samples_breastcancer.Properties.VariableNames)

file = fopen('/Users/s202425/Documents/GitHub/metastasis/results/sampling/breastcancer_vs_lungmetastasis.csv', 'a');
fprintf(file, 'reaction;stats.ks;stats.ranksum;stats.ttest;pVals.ks;pVals.ranksum;pVals.ttest;sampleStats.mean;sampleStats.std;sampleStats.mode;sampleStats.median;sampleStats.skew;sampleStats.kurt');

for pos = 1:length(commonreaction)
   reaction = commonreaction{pos};
   % get sample array from specific array
   samp_breastcancer = samples_breastcancer{:, {reaction}};
   samp_lungmet = samples_lungmet{:, {reaction}};
   % Kolmogorov-Smirnov test, rank-sum test, T-test
   [stats, pVals] = compareTwoSamplesStat(samp_breastcancer,samp_lungmet, {'ks', 'rankSum','tTest'});
   sampleStats = calcSampleStats({samp_breastcancer, samp_lungmet});
   % write to file
   fprintf(file, reaction, ';', mean(stats.ks), ';', mean(stats.ranksum),';', ...
   mean(stats.ttest), ';', mean(pVals.ks), ';', mean(pVals.ranksum), ';', mean(pVals.ttest), ';',...
   sampleStats.mean, ';', sampleStats.std, ';', sampleStats.mode, ';',...
   sampleStats.median, ';',sampleStats.skew, ';',sampleStats.kurt);
end 
fclose(file);
