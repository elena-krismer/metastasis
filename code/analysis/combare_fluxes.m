
% calculate jaccard similarity for flux variabilitz analysis

% brainmetastasis vs. lungtumour
load('/Users/s202425/Documents/GitHub/metastasis/obj/models/lung_metastasis.mat'.mat');
load('/Users/s202425/Documents/GitHub/metastasis/obj/models/brain_metastasis.mat');

% find common reactions in models 
[commonreaction, pos] = intersect(GIMME_group1_brainmetastasis.rxns, GIMME_group2_lungtumour.rxns);
% flux variability analysis 
[minFlux_brainmetastasis, maxFlux_brainmetastasis] = fluxVariability(GIMME_group1_brainmetastasis, 100, "max", commonreactions);
[minFlux_lungtumour, maxFlux_lungtumour] = fluxVariability(GIMME_group2_lungtumour, 100, "max", commonreactions);

% combine min and max fluxes to matrix
min_fluxes =[minFlux_brainmetastasis(:), minFlux_lungtumour(:)];
max_fluxes =[maxFlux_brainmetastasis(:), maxFlux_lungtumour(:)];

% calculate jaccard
J = fvaJaccardIndex(min_fluxes, max_fluxes);
% convert double to cell
J = num2cell(J);
% matrix with reaction name and jaccard 
jaccard_reactions = [commonreactions(:), J(:)]


% breast origin metastais in lung vs. original lungtumour
load('/Users/s202425/Documents/GitHub/metastasis/obj/models/metastasis_lungtissue.mat');
load('/Users/s202425/Documents/GitHub/metastasis/obj/models/lungtumour.mat');
% find common reactions in models
[commonreactions, pos] = intersect(metastasis_lungtissue.rxns, GIMME_group2_lungtumour.rxns);
[minFlux_metastasis, maxFlux_metastasis] = fluxVariability(metastasis_lungtissue, 100, "max", commonreactions);
[minFlux_lungtumour, maxFlux_lungtumour] = fluxVariability(GIMME_group2_lungtumour, 100, "max", commonreactions);
min_fluxes =[minFlux_metastasis(:), minFlux_lungtumour(:)];
max_fluxes =[maxFlux_metastasis(:), maxFlux_lungtumour(:)];
J = fvaJaccardIndex(min_fluxes, max_fluxes);
% convert double to cell
J = num2cell(J);
% matrix with reaction name and jaccard 
jaccard_reactions = [commonreactions(:), J(:)]

% breast origin metastais in lung vs. nonlung
load('/Users/s202425/Documents/GitHub/metastasis/obj/models/metastasis_lungtissue.mat');
load('/Users/s202425/Documents/GitHub/metastasis/obj/models/metastasis_nonlungtissue.mat');
% find common reactions in models
[commonreactions, pos] = intersect(metastasis_lungtissue.rxns, metastasis_nonlungtissue.rxns);
[minFlux_metastasis, maxFlux_metastasis] = fluxVariability(metastasis_lungtissue, 100, "max", commonreactions);
[minFlux_nonlung, maxFlux_nonlung] = fluxVariability(metastasis_nonlungtissue, 100, "max", commonreactions);
min_fluxes =[minFlux_metastasis(:), minFlux_nonlung(:)];
max_fluxes =[maxFlux_metastasis(:), maxFlux_nonlung(:)];
J = fvaJaccardIndex(min_fluxes, max_fluxes);
% convert double to cell
J = num2cell(J);
% matrix with reaction name and jaccard 
jaccard_reactions = [commonreactions(:), J(:)]

% breast origin metastais in lung vs. origin breast cancer
load('/Users/s202425/Documents/GitHub/metastasis/obj/models/MDA_MB_231.mat');
% find common reactions in models
[commonreactions, pos] = intersect(metastasis_lungtissue.rxns, GIMME_lung_GSE2603.rxns);
[minFlux_metastasis, maxFlux_metastasis] = fluxVariability(metastasis_lungtissue, 100, "max", commonreactions);
[minFlux_breastcancer, maxFlux_breastcancer] = fluxVariability(GIMME_lung_GSE2603, 100, "max", commonreactions);
min_fluxes =[minFlux_metastasis(:), minFlux_breastcancer(:)];
max_fluxes =[maxFlux_metastasis(:), maxFlux_breastcancer(:)];
J = fvaJaccardIndex(min_fluxes, max_fluxes);
% convert double to cell
J = num2cell(J);
% matrix with reaction name and jaccard 
jaccard_reactions = [commonreactions(:), J(:)]

% breast origin metastais in nonlungtissue vs. origin breast cancer
load('/Users/s202425/Documents/GitHub/metastasis/obj/models/MDA_MB_231.mat');
% find common reactions in models
[commonreactions, pos] = intersect(metastasis_nonlungtissue.rxns, GIMME_lung_GSE2603.rxns);
[minFlux_metastasis, maxFlux_metastasis] = fluxVariability(metastasis_nonlungtissue, 100, "max", commonreactions);
[minFlux_breastcancer, maxFlux_breastcancer] = fluxVariability(GIMME_lung_GSE2603, 100, "max", commonreactions);
min_fluxes =[minFlux_metastasis(:), minFlux_breastcancer(:)];
max_fluxes =[maxFlux_metastasis(:), maxFlux_breastcancer(:)];
J = fvaJaccardIndex(min_fluxes, max_fluxes);
% convert double to cell
J = num2cell(J);
% matrix with reaction name and jaccard 
jaccard_reactions = [commonreactions(:), J(:)]