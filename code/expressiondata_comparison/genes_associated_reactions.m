high_expressed = {'ENSG00000049541', 'ENSG00000125813', 'ENSG00000165181',...
    'ENSG00000167900', 'ENSG00000130517', 'ENSG00000138326', 'ENSG00000115414',...
    'ENSG00000008988', 'ENSG00000198755', 'ENSG00000184110', 'ENSG00000177954',...
    'ENSG00000177600', 'ENSG00000163468', 'ENSG00000108679', 'ENSG00000099204',...
    'ENSG00000125970', 'ENSG00000026025', 'ENSG00000162889', 'ENSG00000087087', ...
    'ENSG00000117748', 'ENSG00000141556', 'ENSG00000172809', 'ENSG00000004779',...
    'ENSG00000137699', 'ENSG00000133895', 'ENSG00000136717', 'ENSG00000130066',...
    'ENSG00000002726', 'ENSG00000100425', 'ENSG00000012048', 'ENSG00000106089',...
    'ENSG00000197063', 'ENSG00000051180', 'ENSG00000172137', 'ENSG00000179409',...
'ENSG00000189369', 'ENSG00000197299', 'ENSG00000162924', 'ENSG00000010165',...
'ENSG00000128408', 'ENSG00000106682', 'ENSG00000133112', 'ENSG00000112033', ...
'ENSG00000100664', 'ENSG00000114867','ENSG00000104490', 'ENSG00000063046',...
'ENSG00000089157', 'ENSG00000013392', 'ENSG00000166862','ENSG00000125691',...
'ENSG00000103266', 'ENSG00000198736', 'ENSG00000185088', 'ENSG00000070785',...
'ENSG00000185305', 'ENSG00000007047', 'ENSG00000256646', 'ENSG00000122852',...
'ENSG00000130741', 'ENSG00000111752', 'ENSG00000059769', 'ENSG00000144535'}

[Reaclist] = findRxnsActiveWithGenes(brain_metastasis, high_expressed);
% {'HMR_4476'}
[results, ListResults] = findRxnsFromGenes(brain_metastasis, high_expressed);
% {'HMR_4476','m01169c + m01285c + m02039c  <=> m01371c + m01683c ',
% 'Pentose phosphate pathway','','ENSG00000130066 or ENSG00000141504'}

% ENSG00000130066 Diamine acetyltransferase 1 Enzyme which catalyzes the acetylation of polyamines 

% {'HMR_6921','5 m02039m + m02553m + m03103m  -> m02552m + m03102m + 4 m02039i 
% ','Oxidative phosphorylation','','ENSG00000004779 and ENSG00000023228 and 
% ENSG00000065518 and ENSG00000090266 and ENSG00000099795 and ENSG00000109390 
% and ENSG00000110717 and ENSG00000115286 and ENSG00000119013 and ENSG00000119421 
% and ENSG00000125356 and ENSG00000128609 and ENSG00000130414 and ENSG00000131495 
% and ENSG00000136521 and ENSG00000139180 and ENSG00000140990 and ENSG00000145494 
% and ENSG00000147123 and ENSG00000147684 and ENSG00000151366 and ENSG00000158864 
% and ENSG00000160194 and ENSG00000164258 and ENSG00000165264 and ENSG00000166136 
% and ENSG00000167792 and ENSG00000168653 and ENSG00000170906 and ENSG00000174886 
% and ENSG00000178127 and ENSG00000183648 and ENSG00000184752 and ENSG00000184983 
% and ENSG00000186010 and ENSG00000189043 and ENSG00000198695 and ENSG00000198763 
% and ENSG00000198786 and ENSG00000198840 and ENSG00000198886 and ENSG00000198888 
% and ENSG00000212907 and ENSG00000213619'}

% ENSG00000004779 Acyl carrier protein, mitochondrial

% Returns the subsystems associated with the provided genes.
[geneSubSystems,singleList] = findSubSystemsFromGenes(brain_metastasis, high_expressed);
% [' OPadehilnoprstvwxy']
% [' Oadehilnoprstvxy']
% [' Paehnopstwy']

low_expressed = { 'ENSG00000008056', 'ENSG00000105426', 'ENSG00000197959',...
    'ENSG00000108107', 'ENSG00000100353', 'ENSG00000175390', 'ENSG00000109475',...
    'ENSG00000089009', 'ENSG00000100632', 'ENSG00000167658','ENSG00000110958',...
    'ENSG00000171858', 'ENSG00000186468', 'ENSG00000158417' 'ENSG00000204843',...
'ENSG00000100097', 'ENSG00000110955', 'ENSG00000204209', 'ENSG00000185122', ...
'ENSG00000113649', 'ENSG00000063854', 'ENSG00000167768', 'ENSG00000131480',...
'ENSG00000233927', 'ENSG00000227057','ENSG00000138834', 'ENSG00000171431',...
'ENSG00000223501', 'ENSG00000120210', 'ENSG00000241127','ENSG00000124370'}

% Finds all reactions for which the provided genes show sufficient evidence of 
% their absence. (i.e. make the corresponding GPRs always be zero)
[Reaclist] = findRxnsActiveWithGenes(brain_metastasis, low_expressed);
% {'HMR_3857','HMR_8511','HMR_4136','HMR_8604','HMR_4430','HMR_4682',
% 'HMR_3213','HMR_6758','HMR_6789','HMR_7618','HMR_2583','HMR_2584','HMR_8560','HMR_1312','HMR_1394','RE3587N'}
[results, ListResults] = findRxnsFromGenes(brain_metastasis, low_expressed);
% results.gene_ENSG00000110958: 
[geneSubSystems,singleList] = findSubSystemsFromGenes(brain_metastasis, low_expressed);



