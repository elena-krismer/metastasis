% gimme 
% apply GIMME ALGORITHM on expression data

initCobraToolbox(false)
changeCobraSolver ('glpk', 'all');

humanone_model = readCbModel(which('/Users/s202425/Documents/GitHub/Human-GEM/model/Human-GEM.mat'));
options = 'options_GIMME';
load(['options_methods' filesep options]);

pathname = fileparts('/Users/s202425/Documents/GitHub/metastasis/obj/models/');
load('/Users/s202425/Documents/GitHub/metastasis/obj/models/brain_class.mat');

load('/Users/s202425/Documents/GitHub/metastasis/obj/models/lung_metastasis.mat');
load('/Users/s202425/Documents/GitHub/metastasis/obj/models/brain_metastasis.mat');
load('/Users/s202425/Documents/GitHub/metastasis/obj/models/MDA_MB_231.mat');
load('/Users/s202425/Documents/GitHub/metastasis/obj/models/lung_class.mat');
load('/Users/s202425/Documents/GitHub/metastasis/obj/models/breast_class.mat');
load('/Users/s202425/Documents/GitHub/metastasis/obj/models/white_matter_class.mat');
load('/Users/s202425/Documents/GitHub/metastasis/obj/models/gbm_surrounding_tissue_class.mat');
load('/Users/s202425/Documents/GitHub/metastasis/obj/models/glioblastoma_class.mat');
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


[expressionRxns parsedGPR]= mapExpressionToReactions(humanone_model, brain);
% create GIMME model 
options.expressionRxns = expressionRxns;
options.threshold = round(prctile(expressionRxns, 25));
brain_tissue  = createTissueSpecificModel(humanone_model, options);
matfile = fullfile(pathname, 'brain_tissue.mat');
save(matfile, 'brain_tissue');

[expressionRxns parsedGPR]= mapExpressionToReactions(humanone_model, breast);
% create GIMME model 
options.expressionRxns = expressionRxns;
options.threshold = round(prctile(expressionRxns, 25));
breast_tissue  = createTissueSpecificModel(humanone_model, options);
matfile = fullfile(pathname, 'breast_tissue.mat');
save(matfile, 'breast_tissue');

[expressionRxns parsedGPR]= mapExpressionToReactions(humanone_model, lung);
% create GIMME model 
options.expressionRxns = expressionRxns;
options.threshold = round(prctile(expressionRxns, 25));
lung_tissue  = createTissueSpecificModel(humanone_model, options);
matfile = fullfile(pathname, 'lung_tissue.mat');
save(matfile, 'lung_tissue');

[expressionRxns parsedGPR]= mapExpressionToReactions(humanone_model, white_matter);
% create GIMME model 
options.expressionRxns = expressionRxns;
options.threshold = round(prctile(expressionRxns, 25));
white_matter = createTissueSpecificModel(humanone_model, options);
matfile = fullfile(pathname, 'white_matter.mat');
save(matfile, 'white_matter');

[expressionRxns parsedGPR]= mapExpressionToReactions(humanone_model, glioblastoma);
% create GIMME model 
options.expressionRxns = expressionRxns;
options.threshold = round(prctile(expressionRxns, 25));
glioblastoma  = createTissueSpecificModel(humanone_model, options);
matfile = fullfile(pathname, 'glioblastoma.mat');
save(matfile, 'glioblastoma');


[expressionRxns parsedGPR]= mapExpressionToReactions(humanone_model, gbm_surrounding_tissue);
% create GIMME model 
options.expressionRxns = expressionRxns;
options.threshold = round(prctile(expressionRxns, 25));
gbm_surrounding_tissue  = createTissueSpecificModel(humanone_model, options);
matfile = fullfile(pathname, 'gbm_surrounding_tissue.mat');
save(matfile, 'gbm_surrounding_tissue');

% GSE50161
load('/Users/s202425/Documents/GitHub/metastasis/obj/models/glioblastoma_class.mat');
load('/Users/s202425/Documents/GitHub/metastasis/obj/models/medullablastoma_class.mat');
load('/Users/s202425/Documents/GitHub/metastasis/obj/models/pilocyticastrocytoma_class.mat');
load('/Users/s202425/Documents/GitHub/metastasis/obj/models/ependymoma_class.mat');

[expressionRxns parsedGPR]= mapExpressionToReactions(humanone_model, glioblastoma);
% create GIMME model 
options.expressionRxns = expressionRxns;
options.threshold = round(prctile(expressionRxns, 25));
glioblastoma  = createTissueSpecificModel(humanone_model, options);
matfile = fullfile(pathname, 'glioblastoma.mat');
save(matfile, 'glioblastoma');

[expressionRxns parsedGPR]= mapExpressionToReactions(humanone_model, medullablastoma);
% create GIMME model 
options.expressionRxns = expressionRxns;
options.threshold = round(prctile(expressionRxns, 25));
medullablastoma  = createTissueSpecificModel(humanone_model, options);
matfile = fullfile(pathname, 'medullablastoma.mat');
save(matfile, 'medullablastoma');

[expressionRxns parsedGPR]= mapExpressionToReactions(humanone_model,pilocyticastrocytoma);
% create GIMME model 
options.expressionRxns = expressionRxns;
options.threshold = round(prctile(expressionRxns, 25));
pilocyticastrocytoma  = createTissueSpecificModel(humanone_model, options);
matfile = fullfile(pathname, 'pilocyticastrocytoma.mat');
save(matfile, 'pilocyticastrocytoma');

[expressionRxns parsedGPR]= mapExpressionToReactions(humanone_model,ependymoma);
% create GIMME model 
options.expressionRxns = expressionRxns;
options.threshold = round(prctile(expressionRxns, 25));
ependymoma  = createTissueSpecificModel(humanone_model, options);
matfile = fullfile(pathname, 'ependymoma.mat');
save(matfile, 'ependymoma');
