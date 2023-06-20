clear all
%% Spring or fall
% SPRING --------------------------------------------
load 1HertzStats_all_summer2021
SratsT=SratsT(SratsT.MeadianInterval<2,:); % only 1 Hz
%========================================================
%% First time or adding data
%DataAll=[]; % use only first time
%------------ if not the first time ------------------------------
% load DataFor1Hz_updated2
% DataAll=DataAll2;
%--- another option is to make a third file if we adding large amount of data---------
DataAll=[];
load DataFor1Hz_updated1
% upadte StatsDownloaded
StatsDownloaded_old=StatsDownloaded;
StatsDownloaded=SratsT;
StatsDownloaded.FoundData(:,1)=nan;
StatsDownloaded.Ind(:,1)=nan;
for sd=1:height(StatsDownloaded_old)
    indSt=find(SratsT.Individual==StatsDownloaded_old.Individual(sd) &...
        SratsT.TimestartV==StatsDownloaded_old.TimestartV(sd));
    if ~isempty(indSt)
        StatsDownloaded.FoundData(indSt)=StatsDownloaded_old.FoundData(sd);
        StatsDownloaded.Ind(indSt)=StatsDownloaded_old.Ind(sd);
    end
end
%% Use to read more data files
% load DataFor1Hz
% SratsT=StatsDownloaded;
%% Get the files in the direcrory
directory=('../CraneData\1Hz\');
files=dir('../CraneData\1Hz\*.csv');
%% loop over all the files
for ii=1:length(files)
    %% read the csv data for that crane
    tt=readtable([directory,files(ii).name]);
    % check that it's new and wasn't previously added
    In1Hz=find(SratsT.Tag==tt.device_id(1) & ...
        (isnan(StatsDownloaded.FoundData) | StatsDownloaded.FoundData==0)); % find 1Hz pices for this tag
    for h=1:length(In1Hz) % loop for all 1Hz peeces for this tag
        % Only relevant part of data
        SrtT=SratsT.TimestartV(In1Hz(h));
        EndT=SratsT.TimeEndV(In1Hz(h));
        
        disp(['Workin on tag ',num2str(tt.device_id(1)),', starting at ', datestr(datenum(SrtT))])
        
        hj=find(tt.UTC_datetime>SrtT & tt.UTC_datetime<EndT);
        % check if we have this part, maybe it was already done before
        if isempty(hj)
            continue;
        end
        jkkk=diff(hj);
        %tt.device_id(1)==17089 &
        if  ~isempty(find(jkkk>1))
            if length(hj)-(find(jkkk>1))<5 |  length(find(jkkk>1))>1
                Data=tt(tt.UTC_datetime>SrtT & tt.UTC_datetime<EndT,:);
            else
                Data=tt(hj(find(jkkk>1)+1:end),:);
            end
        else
            Data=tt(tt.UTC_datetime>SrtT & tt.UTC_datetime<EndT,:);
        end
        %---------- check if we have time duplicates-----------
        GoingBack=find((diff(datenum(Data.UTC_datetime))*86400)<0);
        gaps=diff(datenum(Data.UTC_datetime))*86400;
        if ~isempty(GoingBack)
            a=[Data.UTC_datetime(GoingBack),Data.UTC_datetime(GoingBack+1)];
            aIndex=[GoingBack,GoingBack+1];
            TempData=Data(1:aIndex(1,1),:);
            for b=1:length(a(:,1))
                indexErazefrom=find(TempData.UTC_datetime==a(b,2));
                if length(indexErazefrom)>1 % if we on ACC data
                    indexErazefrom=indexErazefrom(1);
                end
                if aIndex(b,2)==height(Data)
                    break;
                end
                % if the data is realy duplicated
                if Data.UTC_datetime(aIndex(b,2)+1)==TempData.UTC_datetime(indexErazefrom+1)
                    TempData(indexErazefrom:end,:)=[];
                    if b==length(a(:,1))
                        TempData=[TempData; Data(aIndex(b,2):end,:)];
                    else
                        TempData=[TempData; Data(aIndex(b,2):aIndex(b+1,1),:)];
                    end
                else
                    if b==length(a(:,1))
                        TempData=[TempData; Data(aIndex(b,2)+1:end,:)];
                    else
                        TempData=[TempData; Data(aIndex(b,2)+1:aIndex(b+1,1),:)];
                    end
                    
                end
            end
            Data=TempData;
        end
        %------------------------------------------------------
        Tag=ones(height(Data),1)*tt.device_id(1);
        HZIndex=ones(height(Data),1)*In1Hz(h);
        StatsDownloaded.FoundData(In1Hz(h))=1;
        StatsDownloaded.Ind(In1Hz(h))=In1Hz(h);
        Temp=table(Tag,HZIndex);
        DataAll=[DataAll; [Temp,Data]];
    end
end
% if adding and created another file
DataAll3=DataAll;
StatsDownloaded.FoundData(isnan(StatsDownloaded.FoundData))=0;
StatsDownloaded.Ind(isnan(StatsDownloaded.Ind))=0;
save('DataFor1Hz_updated3','DataAll3','StatsDownloaded')

DataAll1=DataAll(1:3387918,:);
DataAll2=DataAll(3387919:end,:);
save('DataFor1Hz_updated1','DataAll1','StatsDownloaded')
save('DataFor1Hz_updated2','DataAll2','StatsDownloaded')
