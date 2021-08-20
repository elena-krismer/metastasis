initCobraToolbox (0)
changeCobraSolver('gurobi') 

model = readCbModel('sampleMetaOut_EndoRecon1NewBounds.mat');
model_points = model.points; 
Model1_I = readCbModel('Model1_I.mat');
Model2_I = readCbModel('Model2_I.mat');
Model3_I = readCbModel('Model3_I.mat');
Model4_I = readCbModel('Model4_I.mat');
TasksTable=readtable('Supp_Tables_2.xls','sheet','Supp_Table_1_Recon1_Filtered');
Tasks=unique(TasksTable(:,1)); % Hi Elena, this line has to be changed in order to avoid "unique" function sorts the tasks alphabetically. Check "unique" function documentation
ListOfTasks=unique(join(string(TasksTable{:,2:4}),'_'));
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

%Plot
cgo=clustergram(AllTraumaSummaryTable,'Standardize','Row', 'Colormap', 'redbluecmap');
set(cgo,'Linkage','complete','Dendrogram',3, 'ColumnLabels', {'Trauma 1','Trauma 2','Trauma 3','Trauma 4'}, 'ColumnLabelsRotate', 45)
rm = struct('GroupNumber',{180,185,183,184,182},'Annotation',{'Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5'},'Color',{'m','r','g','y','b'});
set(cgo,'RowGroupMarker',rm)

%Plot analysis 1
%Cluster 1: Note: the variable Cluster_1 has been generated manually from %the cgo plot
Pathways_Cluster_1=unique(ListOfTasksAll(cellfun(@str2num,Cluster_1.RowLabels)));
Pathways_Weight_Cluster_1=zeros(length(Pathways_Cluster_1),1)
for i=1:length(Pathways_Cluster_1)
    Pathways_Weight_Cluster_1(i)=sum(ismember(ListOfTasksAll(cellfun(@str2num,Cluster_1.RowLabels)),Pathways_Cluster_1(i)))/length(Cluster_1.RowLabels);
end
%Cluster 2: Note: the variable Cluster_2 has been generated manually from %the cgo plot
Pathways_Cluster_2=unique(ListOfTasksAll(cellfun(@str2num,Cluster_2.RowLabels)));
Pathways_Weight_Cluster_2=zeros(length(Pathways_Cluster_2),1)
for i=1:length(Pathways_Cluster_2)
    Pathways_Weight_Cluster_2(i)=sum(ismember(ListOfTasksAll(cellfun(@str2num,Cluster_2.RowLabels)),Pathways_Cluster_2(i)))/length(Cluster_2.RowLabels);
end
%Cluster 3: Note: the variable Cluster_3 has been generated manually from %the cgo plot
Pathways_Cluster_3=unique(ListOfTasksAll(cellfun(@str2num,Cluster_3.RowLabels)));
Pathways_Weight_Cluster_3=zeros(length(Pathways_Cluster_3),1)
for i=1:length(Pathways_Cluster_3)
    Pathways_Weight_Cluster_3(i)=sum(ismember(ListOfTasksAll(cellfun(@str2num,Cluster_3.RowLabels)),Pathways_Cluster_3(i)))/length(Cluster_3.RowLabels);
end
%Cluster 4: Note: the variable Cluster_4 has been generated manually from %the cgo plot
Pathways_Cluster_4=unique(ListOfTasksAll(cellfun(@str2num,Cluster_4.RowLabels)));
Pathways_Weight_Cluster_4=zeros(length(Pathways_Cluster_4),1)
for i=1:length(Pathways_Cluster_4)
    Pathways_Weight_Cluster_4(i)=sum(ismember(ListOfTasksAll(cellfun(@str2num,Cluster_4.RowLabels)),Pathways_Cluster_4(i)))/length(Cluster_4.RowLabels);
end
%Cluster 5: Note: the variable Cluster_5 has been generated manually from %the cgo plot
Pathways_Cluster_5=unique(ListOfTasksAll(cellfun(@str2num,Cluster_5.RowLabels)));
Pathways_Weight_Cluster_5=zeros(length(Pathways_Cluster_5),1)
for i=1:length(Pathways_Cluster_5)
    Pathways_Weight_Cluster_5(i)=sum(ismember(ListOfTasksAll(cellfun(@str2num,Cluster_5.RowLabels)),Pathways_Cluster_5(i)))/length(Cluster_5.RowLabels);
end

%Plot analysis 2
%Cluster 1: Note: the variable Cluster_1 has been generated manually from %the cgo plot
Pathways_Cluster_1_b=unique(ListOfTasksAll(cellfun(@str2num,Cluster_1.RowLabels),2));
Pathways_Weight_Cluster_1_b=zeros(length(Pathways_Cluster_1_b),1)
for i=1:length(Pathways_Cluster_1_b)
    Pathways_Weight_Cluster_1_b(i)=sum(ismember(ListOfTasksAll(cellfun(@str2num,Cluster_1.RowLabels),2),Pathways_Cluster_1_b(i)))/length(Cluster_1.RowLabels);
end
%Cluster 2: Note: the variable Cluster_2 has been generated manually from %the cgo plot
Pathways_Cluster_2_b=unique(ListOfTasksAll(cellfun(@str2num,Cluster_2.RowLabels),2));
Pathways_Weight_Cluster_2_b=zeros(length(Pathways_Cluster_2_b),1)
for i=1:length(Pathways_Cluster_2_b)
    Pathways_Weight_Cluster_2_b(i)=sum(ismember(ListOfTasksAll(cellfun(@str2num,Cluster_2.RowLabels),2),Pathways_Cluster_2_b(i)))/length(Cluster_2.RowLabels);
end
%Cluster 3: Note: the variable Cluster_3 has been generated manually from %the cgo plot
Pathways_Cluster_3_b=unique(ListOfTasksAll(cellfun(@str2num,Cluster_3.RowLabels),2));
Pathways_Weight_Cluster_3_b=zeros(length(Pathways_Cluster_3_b),1)
for i=1:length(Pathways_Cluster_3_b)
    Pathways_Weight_Cluster_3_b(i)=sum(ismember(ListOfTasksAll(cellfun(@str2num,Cluster_3.RowLabels),2),Pathways_Cluster_3_b(i)))/length(Cluster_3.RowLabels);
end
%Cluster 4: Note: the variable Cluster_4 has been generated manually from %the cgo plot
Pathways_Cluster_4_b=unique(ListOfTasksAll(cellfun(@str2num,Cluster_4.RowLabels),2));
Pathways_Weight_Cluster_4_b=zeros(length(Pathways_Cluster_4_b),1)
for i=1:length(Pathways_Cluster_4_b)
    Pathways_Weight_Cluster_4_b(i)=sum(ismember(ListOfTasksAll(cellfun(@str2num,Cluster_4.RowLabels),2),Pathways_Cluster_4_b(i)))/length(Cluster_4.RowLabels);
end
%Cluster 5: Note: the variable Cluster_5 has been generated manually from %the cgo plot
Pathways_Cluster_5_b=unique(ListOfTasksAll(cellfun(@str2num,Cluster_5.RowLabels),2));
Pathways_Weight_Cluster_5_b=zeros(length(Pathways_Cluster_5_b),1)
for i=1:length(Pathways_Cluster_5_b)
    Pathways_Weight_Cluster_5_b(i)=sum(ismember(ListOfTasksAll(cellfun(@str2num,Cluster_5.RowLabels),2),Pathways_Cluster_5_b(i)))/length(Cluster_5.RowLabels);
end


ListUniqueTask=unique(ListOfTasksAll(:,2));
SummaryTableTaskLevel_2=zeros(length(ListUniqueTask),5);
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
    if sum(ismember(Pathways_Cluster_5_b,ListUniqueTask(i)))~=0  
        SummaryTableTaskLevel_2(i,5)=Pathways_Weight_Cluster_5_b(find(ismember(Pathways_Cluster_5_b,ListUniqueTask(i)))); % In cluster5
    else
        SummaryTableTaskLevel_2(i,5)=0;
    end 
end
