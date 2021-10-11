solverName = 'glpk';
solverType = 'LP';
changeCobraSolver(solverName, solverType);

load('/Users/s202425/Documents/GitHub/metastasis/obj/models/lung_metastasis.mat')
load('/Users/s202425/Documents/GitHub/metastasis/obj/models/brain_metastasis.mat');

%protons m02039c

[inform, m, model_lung] = checkStoichiometricConsistency(lung_metastasis)
[inform, m, model_breast] = checkStoichiometricConsistency(MDA_MB_231)
[inform, m, model_brain] = checkStoichiometricConsistency(brain_metastasis)

[leakMetBool_lung, leakRxnBool_lung, siphonMetBool_lung, siphonRxnBool_lung, leakY_lung, siphonY_lung, statp_lung, statn_lung] = findMassLeaksAndSiphons(model_lung)
% leakMetBool: 821 x 1
% leakRxnBool: 740 x 0/1
% siphonMetBool 821 x1
% siphonRxnBool_lung 740 x 0/1
% leakY_lung: (409—414,1) 0.2634
% siphonY_lung: (409—414,1) 0.2634
% m02039c, m02039m, m02039n, m02039p, m02039r, m02039s



[leakMetBool_brain, leakRxnBool_brain, siphonMetBool_brain, siphonRxnBool_brain, leakY_brain, siphonY_brain, statp_brain, statn_brain] = findMassLeaksAndSiphons(model_brain)
% leakMetBool: 813 x 1
% leakRxnBool: 734 x 0/1
% siphonMetBool 813 x1
% siphonRxnBool 734 x 0/1
% leakY: (402—407,1) 0.2634
% siphonY: (402—407,1) 0.2634
% m02039c, m02039m, m02039n, m02039p, m02039r, m02039s


[leakMetBool_breast, leakRxnBool_breast, siphonMetBool_breast, siphonRxnBool_breast, leakY_breast, siphonY_breast, statp_breast statn_breast] = findMassLeaksAndSiphons(model_breast)
% leakMetBool: 809 x 1
% leakRxnBool: 729 x 0/1
% siphonMetBool 809 x1
% siphonRxnBool 729 x 0/1
% leakY: (401—406,1) 0.2634
% siphonY: (401—406,1) 0.2634
% m02039c, m02039m, m02039n, m02039p, m02039r, m02039s