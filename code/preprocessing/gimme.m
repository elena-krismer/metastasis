% gimme 
% apply GIMME ALGORITHM on expression data

initCobraToolbox(false)
changeCobraSolver ('glpk', 'all');

humanone_model = readCbModel(which('/Users/s202425/Documents/GitHub/breast_metastasis/obj/models/Human-GEM.mat'));
options = 'options_GIMME';
load(['options_methods' filesep options]);

pathname = fileparts('/Users/s202425/Documents/GitHub/metastasis/obj/models/');


load('/Users/s202425/Documents/GitHub/metastasis/obj/models/lung_metastasis.mat');
load('/Users/s202425/Documents/GitHub/metastasis/obj/models/brain_metastasis.mat');
load('/Users/s202425/Documents/GitHub/metastasis/obj/models/MDA_MB_231.mat');

% map expression to reaction
% create GIMME model 
[expressionRxns parsedGPR]= mapExpressionToReactions(humanone_model, lung_metastasis);
options.expressionRxns = expressionRxns;
options.threshold = round(prctile(expressionRxns, 25));
lung_metastasis = createTissueSpecificModel(humanone_model, options);
matfile = fullfile(pathname, 'lung_metastasis.mat');
save(matfile, 'lung_metastasis');


% map expression to reaction
[expressionRxns parsedGPR]= mapExpressionToReactions(humanone_model, brain_metastasis);
% create GIMME model 
options.expressionRxns = expressionRxns;
options.threshold = round(prctile(expressionRxns, 25));
brain_metastasis  = createTissueSpecificModel(humanone_model, options);
matfile = fullfile(pathname, 'brain_metastasis.mat');
save(matfile, 'brain_metastasis');

[expressionRxns parsedGPR]= mapExpressionToReactions(humanone_model, MDA_MB_231);
% create GIMME model 
options.expressionRxns = expressionRxns;
options.threshold = round(prctile(expressionRxns, 25));
MDA_MB_231  = createTissueSpecificModel(humanone_model, options);
matfile = fullfile(pathname, 'MDA_MB_231.mat');
save(matfile, 'MDA_MB_231');
