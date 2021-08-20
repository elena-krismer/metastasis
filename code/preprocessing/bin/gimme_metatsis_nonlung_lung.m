% apply GIMME ALGORITHM on expression data

initCobraToolbox(false)
changeCobraSolver ('glpk', 'all');
% load humanone model
humanone_model = readCbModel(which('/Users/s202425/Documents/GitHub/breast_metastasis/obj/models/Human-GEM.mat'));
options = 'options_GIMME';
load(['options_methods' filesep options]);
tissuespecificmodel = createTissueSpecificModel(humanone_model, options);

% load expressiondata
load('/Users/s202425/Documents/GitHub/breast_metastasis/obj/models/model_metastasis_lungtissue.mat');
load('/Users/s202425/Documents/GitHub/breast_metastasis/obj/models/model_metastasis_nonlungtissue.mat');

% map expression to reaction
%expressiondata_GSE11078.genes = expressiondata_GSE11078.gene;
[expressionRxns parsedGPR]= mapExpressionToReactions(humanone_model, model_metastasis_lungtissue);
% create GIMME model 
threshold=round(prctile(expressionRxns, 25));
metastasis_lungtissue = GIMME(tissuespecificmodel, expressionRxns, threshold);

pathname = fileparts('/Users/s202425/Documents/GitHub/metastasis/obj/models/');
matfile = fullfile(pathname, 'metastasis_lungtissue.mat');
save(matfile, 'metastasis_lungtissue');

% map expression to reaction
[expressionRxns parsedGPR]= mapExpressionToReactions(humanone_model, model_metastasis_nonlungtissue);
% create GIMME model 
threshold=round(prctile(expressionRxns, 25));
metastasis_nonlungtissue = GIMME(tissuespecificmodel, expressionRxns, threshold);
matfile = fullfile(pathname, 'metastasis_nonlungtissue.mat');
save(matfile, 'metastasis_nonlungtissue')


