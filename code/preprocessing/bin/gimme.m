% gimme 
% load barcelona data into gimme

humanone_model = readCbModel(which('/Users/s202425/Documents/GitHub/breast_metastasis/obj/models/Human-GEM.mat'));
options = 'options_GIMME';
load(['options_methods' filesep options]);
tissuespecificmodel = createTissueSpecificModel(humanone_model, options);

pathname = fileparts('/Users/s202425/Documents/GitHub/metastasis/obj/models/');
threshold=round(prctile(expressionRxns, 25));
load('/Users/s202425/Documents/GitHub/metastasis/obj/models/primary_breast_cancer.mat');
load('/Users/s202425/Documents/GitHub/metastasis/obj/models/lung_metastasis.mat');

threshold=round(prctile(expressionRxns, 25));
primary_breast_cancer = GIMME(tissuespecificmodel, expressionRxns, threshold);
matfile = fullfile(pathname, 'primary_breast_cancer.mat');
save(matfile, 'primary_breast_cancer');

% map expression to reaction
[expressionRxns parsedGPR]= mapExpressionToReactions(humanone_model, lung_metastasis);
% create GIMME model 
threshold=round(prctile(expressionRxns, 25));
lung_metastasis = GIMME(tissuespecificmodel, expressionRxns, threshold);
matfile = fullfile(pathname, 'lung_metastasis.mat');
save(matfile, 'lung_metastasis');

 
load('/Users/s202425/Documents/GitHub/metastasis/obj/models/lungtissue.mat');
% map expression to reaction
[expressionRxns parsedGPR]= mapExpressionToReactions(humanone_model, lungtissue);
% create GIMME model 
threshold=round(prctile(expressionRxns, 25));
lungtissue = GIMME(tissuespecificmodel, expressionRxns, threshold);
matfile = fullfile(pathname, 'lungtissue.mat');
save(matfile, 'lungtissue');