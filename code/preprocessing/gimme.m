% gimme 
% apply GIMME ALGORITHM on expression data

initCobraToolbox(false)
changeCobraSolver ('glpk', 'all');

humanone_model = readCbModel(which('/Users/s202425/Documents/GitHub/breast_metastasis/obj/models/Human-GEM.mat'));
options = 'options_GIMME';
load(['options_methods' filesep options]);
tissuespecificmodel = createTissueSpecificModel(humanone_model, options);

pathname = fileparts('/Users/s202425/Documents/GitHub/metastasis/obj/models/');

load('/Users/s202425/Documents/GitHub/metastasis/obj/models/XXXX_brain_metastasis.mat');
load('/Users/s202425/Documents/GitHub/metastasis/obj/models/lung_metastasis.mat');

load('/Users/s202425/Documents/GitHub/metastasis/obj/models/brain_metastasis.mat');
% map expression to reaction
[expressionRxns parsedGPR]= mapExpressionToReactions(humanone_model, brain_metastasis);
% create GIMME model 
threshold=round(prctile(expressionRxns, 25));
XXXX_brain_metastasis = GIMME(tissuespecificmodel, expressionRxns, 0.5);
matfile = fullfile(pathname, 'XXXX_brain_metastasis.mat');
save(matfile, 'XXXX_brain_metastasis');


% map expression to reaction
[expressionRxns parsedGPR]= mapExpressionToReactions(humanone_model, lung_metastasis);
% create GIMME model 
threshold=round(prctile(expressionRxns, 25));
lung_metastasis = GIMME(tissuespecificmodel, expressionRxns, threshold);
matfile = fullfile(pathname, 'lung_metastasis.mat');
save(matfile, 'lung_metastasis');

 
load('/Users/s202425/Documents/GitHub/metastasis/obj/models/brain_metastasis.mat');
% map expression to reaction
[expressionRxns parsedGPR]= mapExpressionToReactions(humanone_model, brain_metastasis);
% create GIMME model 
threshold=round(prctile(expressionRxns, 25));
brain_metastasis = GIMME(tissuespecificmodel, expressionRxns, threshold);
matfile = fullfile(pathname, 'brain_metastasis.mat');
save(matfile, 'brain_metastasis');