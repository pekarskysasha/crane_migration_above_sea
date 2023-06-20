clear all
update=0;
addpath('ThermalAnalysis')
addpath('ThermalAnalysis/Working folder/FinalTracks')
%update=1;
if update==0
    load StatsTCleaned % meta data of the clened parts
    Indiv=unique(StatsTCleaned.Individual);
else
    load StatsTCleaned_2update % meta data of the parts to update
    Indiv=unique(StatsTCleaned_2add.Individual);
    StatsTCleaned=StatsTCleaned_2add;
end
%% the directory with the data per crane
if update==0
    directory='ThermalAnalysis/Working folder/FinalTracks';
else
    directory='ThermalAnalysis/Working folder/FinalTracks/Updated';
end
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


FlapAll=[];
SG_all=[];
TheramlStatsAll=[];
for i=1:length(Indiv)
    Animal_ID=Indiv(i);
    indArchaive=find(IndevidualsExist==Animal_ID);
    if isempty(indArchaive)
        continue;
    end
    load([directory,'\1HzGPS_and_ACC_annotated_',num2str(IndevidualsExist(indArchaive)),'.mat']);
    %TableGPSAllNew=eval(['TableGPSAllNew_',mat2str(Animal_ID)]);
    AllInfo=StatsTCleaned(StatsTCleaned.Individual==Indiv(i),:);
    day_lightAll=[];
     for h=1:height(AllInfo)
        disp(['working on individual: ',num2str(Animal_ID),' section starts on ',char(AllInfo.TimeStartV(h))])
        TableGPSCont=TableGPSAllNew(TableGPSAllNew.UniqueSectionCounter==AllInfo.UniqueSectionCounter(h),:);
        
        [SunRiseSet]=suncycle(TableGPSCont.Latitude(1),TableGPSCont.Longitude(1),datenum(TableGPSCont.TimeCont(1)));
        SunRiseSet2use=SunRiseSet/24; %sunrise and sunset in UTC
        sunriseThisDay=SunRiseSet2use(1)+floor(datenum(TableGPSCont.TimeCont(1)));
        sunsetThisDay=SunRiseSet2use(2)+floor(datenum(TableGPSCont.TimeCont(1)));
        day_light=datenum(TableGPSCont.TimeCont)>=sunriseThisDay &...
            datenum(TableGPSCont.TimeCont)< sunsetThisDay;
        day_lightAll=[day_lightAll; day_light];
     end
     TableGPSAllNew.day_light=day_lightAll;
     save([directory,'\1HzGPS_and_ACC_annotated_',num2str(IndevidualsExist(indArchaive))],'TableGPSAllNew')
end