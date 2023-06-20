clear all
addpath('ThermalAnalysis')
addpath('ThermalAnalysis/Working folder/FinalTracks')
addpath('ThermalAnalysis/Working folder/FinalTracks/Updated')
 
load StatsTCleaned_2update % meta data of the parts to update
Indiv=unique(StatsTCleaned_2add.Individual);
%% the directory with the data per crane
directory='ThermalAnalysis/Working folder/FinalTracks/Updated';
FilesIn=struct2table(dir(directory));
IndevidualsExistTemp=FilesIn.name(FilesIn.bytes>0);
IndevidualsExist=[];
for ind=1:length(IndevidualsExistTemp)
    dd=char(IndevidualsExistTemp(ind));
    if strcmp(dd(end-2:end),'mat')==0
        continue;
    end
    IndevidualsExist(ind,1)=str2double(dd(26:end-4));
end

directory_main='ThermalAnalysis/Working folder/FinalTracks';
FilesIn=struct2table(dir(directory_main));
IndevidualsExistTemp=FilesIn.name(FilesIn.bytes>0);
IndevidualsExist_old=[];
for ind=1:length(IndevidualsExistTemp)
    dd=char(IndevidualsExistTemp(ind));
    if strcmp(dd(end-2:end),'mat')==0
        continue;
    end
    IndevidualsExist_old(ind,1)=str2double(dd(26:end-4));
end

for i=1:length(Indiv)
    Animal_ID=Indiv(i);
    indArchaive=find(IndevidualsExist==Animal_ID);
    if isempty(indArchaive)
        continue;
    end
    % load the data for the update
    load([directory,'\1HzGPS_and_ACC_annotated_',num2str(IndevidualsExist(indArchaive)),'.mat']);
    TableGPSAllNew_2add=[TableGPSAllNew(:,1:27),TableGPSAllNew(:,30),TableGPSAllNew(:,29),TableGPSAllNew(:,31:end)];   
    % does this indevidual already present in the database
    indArchaive_old=find(IndevidualsExist_old==Animal_ID);
    if ~isempty(IndevidualsExist_old)
        load([directory_main,'\1HzGPS_and_ACC_annotated_',num2str(IndevidualsExist_old(indArchaive_old)),'.mat']);
        TableGPSAllNew=[TableGPSAllNew; TableGPSAllNew_2add];
        TableGPSAllNew=sortrows(TableGPSAllNew,{'TimeCont'});
        save([directory_main,'\1HzGPS_and_ACC_annotated_',num2str(IndevidualsExist_old(indArchaive_old))],'TableGPSAllNew')
    else
        TableGPSAllNew=TableGPSAllNew_2add;
        save([directory_main,'\1HzGPS_and_ACC_annotated_',num2str(IndevidualsExist_old(indArchaive_old))],'TableGPSAllNew')
    end
end
load StatsTCleaned
StatsTCleaned=[StatsTCleaned; StatsTCleaned_2add];
StatsTCleaned=sortrows(StatsTCleaned,{'Individual','UniqueSectionCounter'});
save('StatsTCleaned','StatsTCleaned')