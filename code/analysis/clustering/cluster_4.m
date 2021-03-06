initCobraToolbox (0) 
changeCobraSolver ('glpk', 'all');

Model1_I = readCbModel('/Users/s202425/Documents/GitHub/metastasis/obj/models/lung_metastasis.mat');
Model2_I = readCbModel('/Users/s202425/Documents/GitHub/metastasis/obj/models/brain_metastasis.mat');
Model3_I = readCbModel('/Users/s202425/Documents/GitHub/metastasis/obj/models/MDA_MB_231.mat');
Model4_I = readCbModel('/Users/s202425/Documents/GitHub/metastasis/obj/models/brain_tissue.mat');

TasksTable=readtable('/Users/s202425/Documents/GitHub/metastasis/data/subsystems_tables/Supp_Tables_2_ID_translation_translated.csv');
Tasks=unique(TasksTable(:,1), 'stable'); 
ListOfTasks=unique(join(string(TasksTable{:,2:4}),'_'), 'stable');
ListOfTasksSplit=split(ListOfTasks,'_');
AllTraumaSummaryTable=zeros(height(Tasks),4);

for i=1:height(Tasks)
    %Ith Function
    IthIndex=find(ismember(TasksTable(:,1),Tasks(i,1)));  
    IthSubstrate=table2cell(TasksTable(IthIndex,5));
    IthProduct=table2cell(TasksTable(IthIndex,8));
    %Testing trauma groups
    [Flux_1, FBAsolution_1,~]=testPathway(Model1_I,IthSubstrate,IthProduct);
    [Flux_2, FBAsolution_2,~]=testPathway(Model2_I,IthSubstrate,IthProduct);
    [Flux_3, FBAsolution_3,~]=testPathway(Model3_I,IthSubstrate,IthProduct);   
    [Flux_4, FBAsolution_4,~]=testPathway(Model4_I,IthSubstrate,IthProduct);
    AllTraumaSummaryTable(i,1)=Flux_1;
    AllTraumaSummaryTable(i,2)=Flux_2;
    AllTraumaSummaryTable(i,3)=Flux_3;
    AllTraumaSummaryTable(i,4)=Flux_4;
end

cgo =clustergram(AllTraumaSummaryTable,'Standardize','Row', 'Colormap', 'redbluecmap');
set(cgo,'Linkage','complete','Dendrogram',4, 'ColumnLabels', {'Lung metastasis','Brain metastasis','Breast Cancer', 'Brain Tissue'}, 'ColumnLabelsRotate', 45)
set(cgo, 'RowLabels', cellstr(unique(ListOfTasks(cellfun(@str2num,cgo.RowLabels)), 'stable')))

energy_metabolism_table = AllTraumaSummaryTable(1:8,:);
nucleotide_table = AllTraumaSummaryTable(9:25,:);
carbohydrate_table = AllTraumaSummaryTable(26:60,:);
aminoacid_table = AllTraumaSummaryTable(61:133,:);
lipid_table = AllTraumaSummaryTable(134:171,:);
vitamine_table = AllTraumaSummaryTable(172:178,:);
glycan_table = AllTraumaSummaryTable(179:190);

set(energy_metabolism,'Linkage','complete','Dendrogram',3, 'ColumnLabels', {'Lung metastasis','Brain metastasis','Breast Cancer'}, 'ColumnLabelsRotate', 45);

energy_metabolism=clustergram(energy_metabolism_table,'Standardize','Row', 'Colormap', 'redbluecmap');
nucleotide_metabolism=clustergram(nucleotide_table,'Standardize','Row', 'Colormap', 'redbluecmap');
carbohydrate_metabolism=clustergram(carbohydrate_table,'Standardize','Row', 'Colormap', 'redbluecmap');
aminoacid_metabolism=clustergram(aminoacid_table,'Standardize','Row', 'Colormap', 'redbluecmap');
lipid_metabolism=clustergram(lipid_table,'Standardize','Row', 'Colormap', 'redbluecmap');
vitamine_metabolism=clustergram(vitamine_table,'Standardize','Row', 'Colormap', 'redbluecmap');
glycan_metabolism=clustergram(glycan_table,'Standardize','Row', 'Colormap', 'redbluecmap');



nucleotide_tasks = readtable('GitHub/metastasis/data/pls_da_analysis/nucleotide_tasks.csv');
nucleotide_tasks=join(string(nucleotide_tasks{:,:}),' ');
set(nucleotide_metabolism, 'RowLabels', cellstr(unique(nucleotide_tasks(cellfun(@str2num,nucleotide_metabolism.RowLabels)), 'stable')));
set(nucleotide_metabolism,'Linkage','complete','Dendrogram',4, 'ColumnLabels', {'Lung metastasis','Brain metastasis','Breast Cancer', 'Brain Tissue'}, 'ColumnLabelsRotate', 45);

carbohydrate_tasks = readtable('GitHub/metastasis/data/pls_da_analysis/carbohydrate_tasks.csv');
carbohydrate_tasks=join(string(carbohydrate_tasks{:,:}),' ');
set(carbohydrate_metabolism, 'RowLabels', cellstr(unique(ccarbohydrate_tasks(cellfun(@str2num,carbohydrate_metabolism.RowLabels)), 'stable')));
set(carbohydrate_metabolism, 'RowLabels',cellstr(carbohydrate_tasks(cellfun(@str2num,carbohydrate_metabolism.RowLabels))));
set(carbohydrate_metabolism,'Linkage','complete','Dendrogram',4, 'ColumnLabels', {'Lung metastasis','Brain metastasis','Breast Cancer', 'Brain Tissue'}, 'ColumnLabelsRotate', 45);

aminoacid_tasks = readtable('GitHub/metastasis/data/pls_da_analysis/aminoacid_tasks.csv');
aminoacid_tasks=join(string(aminoacid_tasks{:,:}),'_');
set(aminoacid_metabolism, 'RowLabels', cellstr(unique(aminoacid_tasks(cellfun(@str2num,aminoacid_metabolism.RowLabels)), 'stable')));
set(aminoacid_metabolism,'Linkage','complete','Dendrogram',3, 'ColumnLabels', {'Lung metastasis','Brain metastasis','Breast Cancer'}, 'ColumnLabelsRotate', 45);

lipid_tasks = readtable('GitHub/metastasis/data/pls_da_analysis/lipid_tasks.csv');
lipid_tasks=join(string(lipid_tasks{:,:}),' ');
set(lipid_metabolism, 'RowLabels', cellstr(unique(lipid_tasks(cellfun(@str2num,lipid_metabolism.RowLabels)), 'stable')));
set(lipid_metabolism,'Linkage','complete','Dendrogram',3, 'ColumnLabels', {'Lung metastasis','Brain metastasis','Breast Cancer'}, 'ColumnLabelsRotate', 45);

vitamine_tasks = readtable('GitHub/metastasis/data/pls_da_analysis/vitamine_tasks.csv');
vitamine_tasks=join(string(vitamine_tasks{:,:}),' ');
set(vitamine_metabolism, 'RowLabels', cellstr(unique(vitamine_tasks(cellfun(@str2num,vitamine_metabolism.RowLabels)), 'stable')));
set(vitamine_metabolism,'Linkage','complete','Dendrogram',3, 'ColumnLabels', {'Lung metastasis','Brain metastasis','Breast Cancer'}, 'ColumnLabelsRotate', 45);

glycan_tasks = readtable('GitHub/metastasis/data/pls_da_analysis/glycan_tasks.csv');
glycan_tasks=join(string(glycan_tasks{:,:}),'_');
set(glycan_metabolism, 'RowLabels', cellstr(unique(glycan_tasks(cellfun(@str2num,glycan_metabolism.RowLabels)), 'stable')));
set(nglycan_metabolism,'Linkage','complete','Dendrogram',3, 'ColumnLabels', {'Lung metastasis','Brain metastasis','Breast Cancer'}, 'ColumnLabelsRotate', 45);


%Plot
cgo_1=clustergram(AllTraumaSummaryTable,'Standardize','Row', 'Colormap', 'redbluecmap');
set(cgocgo_1,'Linkage','complete','Dendrogram',3, 'ColumnLabels', {'Lung metastasis','Brain metastasis','Breast Cancer'}, 'ColumnLabelsRotate', 45)
% find nodes where different clusters are divided 
%rm = struct('GroupNumber',{180,185,183,184,182},'Annotation',{'Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5'},'Color',{'m','r','g','y','b'});
%set(cgo,'RowGroupMarker',rm)
rm = struct('GroupNumber',{180,185,186,182},'Annotation',{'Cluster 1','Cluster 2','Cluster 3','Cluster 4'},'Color',{'m','r','g','y'});
set(cgo,'RowGroupMarker',rm)
set(cgo, 'RowLabels', cellstr(unique(ListOfTasks(cellfun(@str2num,cgo.RowLabels)), 'stable')))
%Plot analysis 1
%Cluster 1: Note: the variable Cluster_1 has been generated manually from %the cgo plot
Pathways_Cluster_1=unique(ListOfTasks(cellfun(@str2num,Cluster_1.RowLabels)), 'stable');

Pathways_Weight_Cluster_1=zeros(length(Pathways_Cluster_1),1)
for i=1:length(Pathways_Cluster_1)
    Pathways_Weight_Cluster_1(i)=sum(ismember(ListOfTasks(cellfun(@str2num,Cluster_1.RowLabels)),Pathways_Cluster_1(i)))/length(Cluster_1.RowLabels);
end
%Cluster 2: Note: the variable Cluster_2 has been generated manually from %the cgo plot
Pathways_Cluster_2=unique(ListOfTasks(cellfun(@str2num,Cluster_2.RowLabels)), 'stable');
Cluster_2.RowLabels = Pathways_Cluster_2
Pathways_Weight_Cluster_2=zeros(length(Pathways_Cluster_2),1)
for i=1:length(Pathways_Cluster_2)
    Pathways_Weight_Cluster_2(i)=sum(ismember(ListOfTasks(cellfun(@str2num,Cluster_2.RowLabels)),Pathways_Cluster_2(i)))/length(Cluster_2.RowLabels);
end
%Cluster 3: Note: the variable Cluster_3 has been generated manually from %the cgo plot
Pathways_Cluster_3=unique(ListOfTasks(cellfun(@str2num,Cluster_3.RowLabels)), 'stable');
Pathways_Weight_Cluster_3=zeros(length(Pathways_Cluster_3),1)
for i=1:length(Pathways_Cluster_3)
    Pathways_Weight_Cluster_3(i)=sum(ismember(ListOfTasks(cellfun(@str2num,Cluster_3.RowLabels)),Pathways_Cluster_3(i)))/length(Cluster_3.RowLabels);
end
%Cluster 4: Note: the variable Cluster_4 has been generated manually from %the cgo plot
Pathways_Cluster_4=unique(ListOfTasks(cellfun(@str2num,Cluster_4.RowLabels)), 'stable');
Pathways_Weight_Cluster_4=zeros(length(Pathways_Cluster_4),1)
for i=1:length(Pathways_Cluster_4)
    Pathways_Weight_Cluster_4(i)=sum(ismember(ListOfTasks(cellfun(@str2num,Cluster_4.RowLabels)),Pathways_Cluster_4(i)))/length(Cluster_4.RowLabels);
end
%Cluster 5: Note: the variable Cluster_5 has been generated manually from %the cgo plot
Pathways_Cluster_5=unique(ListOfTasksAll(cellfun(@str2num,Cluster_5.RowLabels)), 'stable');
Pathways_Weight_Cluster_5=zeros(length(Pathways_Cluster_5),1)
for i=1:length(Pathways_Cluster_5)
    Pathways_Weight_Cluster_5(i)=sum(ismember(ListOfTasksAll(cellfun(@str2num,Cluster_5.RowLabels)),Pathways_Cluster_5(i)))/length(Cluster_5.RowLabels);
end

%Plot analysis 2
%Cluster 1: Note: the variable Cluster_1 has been generated manually from %the cgo plot
Pathways_Cluster_1_b=unique(ListOfTasksSplit(cellfun(@str2num,Cluster_1.RowLabels),2), 'stable');
Pathways_Weight_Cluster_1_b=zeros(length(Pathways_Cluster_1_b),1)
for i=1:length(Pathways_Cluster_1_b)
    Pathways_Weight_Cluster_1_b(i)=sum(ismember(ListOfTasksSplit(cellfun(@str2num,Cluster_1.RowLabels),2),Pathways_Cluster_1_b(i)))/length(Cluster_1.RowLabels);
end
%Cluster 2: Note: the variable Cluster_2 has been generated manually from %the cgo plot
Pathways_Cluster_2_b=unique(ListOfTasksSplit(cellfun(@str2num,Cluster_2.RowLabels),2), 'stable');
Pathways_Weight_Cluster_2_b=zeros(length(Pathways_Cluster_2_b),1)
for i=1:length(Pathways_Cluster_2_b)
    Pathways_Weight_Cluster_2_b(i)=sum(ismember(ListOfTasksSplit(cellfun(@str2num,Cluster_2.RowLabels),2),Pathways_Cluster_2_b(i)))/length(Cluster_2.RowLabels);
end
%Cluster 3: Note: the variable Cluster_3 has been generated manually from %the cgo plot
Pathways_Cluster_3_b=unique(ListOfTasksSplit(cellfun(@str2num,Cluster_3.RowLabels),2),'stable');
Pathways_Weight_Cluster_3_b=zeros(length(Pathways_Cluster_3_b),1)
for i=1:length(Pathways_Cluster_3_b)
    Pathways_Weight_Cluster_3_b(i)=sum(ismember(ListOfTasksSplit(cellfun(@str2num,Cluster_3.RowLabels),2),Pathways_Cluster_3_b(i)))/length(Cluster_3.RowLabels);
end
%Cluster 4: Note: the variable Cluster_4 has been generated manually from %the cgo plot
Pathways_Cluster_4_b=unique(ListOfTasksSplit(cellfun(@str2num,Cluster_4.RowLabels),2), 'stable');
Pathways_Weight_Cluster_4_b=zeros(length(Pathways_Cluster_4_b),1)
for i=1:length(Pathways_Cluster_4_b)
    Pathways_Weight_Cluster_4_b(i)=sum(ismember(ListOfTasksSplit(cellfun(@str2num,Cluster_4.RowLabels),2),Pathways_Cluster_4_b(i)))/length(Cluster_4.RowLabels);
end
%Cluster 5: Note: the variable Cluster_5 has been generated manually from %the cgo plot
Pathways_Cluster_5_b=unique(ListOfTasksAll(cellfun(@str2num,Cluster_5.RowLabels),2), 'stable');
Pathways_Weight_Cluster_5_b=zeros(length(Pathways_Cluster_5_b),1)
for i=1:length(Pathways_Cluster_5_b)
    Pathways_Weight_Cluster_5_b(i)=sum(ismember(ListOfTasksAll(cellfun(@str2num,Cluster_5.RowLabels),2),Pathways_Cluster_5_b(i)))/length(Cluster_5.RowLabels);
end


% ListUniqueTask=unique(ListOfTasksAll(:,2));
% check documentation not order 
ListUniqueTask=unique(ListOfTasks(:,1), 'stable');
SummaryTableTaskLevel_2=zeros(length(ListUniqueTask),4);
for i=1:length(ListUniqueTask)
    if sum(ismember(Pathways_Cluster_1_b,ListUniqueTask(i)))~=0
        SummaryTableTaskLevel_2(i,1)=Pathways_Weight_Cluster_1_b(find(ismember(Pathways_Cluster_1_b,ListUniqueTask(i)))); % In cluster1
    else
        SummaryTableTaskLevel_2(i,1)=0;
    end
    if sum(ismember(Pathways_Cluster_2_b,ListUniqueTask(i)))~=0
        SummaryTableTaskLevel_2(i,2)=Pathways_Weight_Cluster_2_b(find(ismember(Pathways_Cluster_2_b,ListUniqueTask(i)))); % In cluster2
    else
        SummaryTableTaskLevel_2(i,2)=0;
    end
    if sum(ismember(Pathways_Cluster_3_b,ListUniqueTask(i)))~=0
        SummaryTableTaskLevel_2(i,3)=Pathways_Weight_Cluster_3_b(find(ismember(Pathways_Cluster_3_b,ListUniqueTask(i)))); % In cluster3
    else
        SummaryTableTaskLevel_2(i,3)=0;
    end   
    if sum(ismember(Pathways_Cluster_4_b,ListUniqueTask(i)))~=0
        SummaryTableTaskLevel_2(i,4)=Pathways_Weight_Cluster_4_b(find(ismember(Pathways_Cluster_4_b,ListUniqueTask(i)))); % In cluster4
    else
        SummaryTableTaskLevel_2(i,4)=0;
    end           
    %if sum(ismember(Pathways_Cluster_5_b,ListUniqueTask(i)))~=0  
       % SummaryTableTaskLevel_2(i,5)=Pathways_Weight_Cluster_5_b(find(ismember(Pathways_Cluster_5_b,ListUniqueTask(i)))); % In cluster5
   % else
   %     SummaryTableTaskLevel_2(i,5)=0;
   % end 
end
