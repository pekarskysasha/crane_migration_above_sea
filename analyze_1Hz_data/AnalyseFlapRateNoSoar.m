function [Answer]=AnalyseFlapRateNoSoar(ThermalTable,Lat,Lon,T,ElevationEboveG,FR,OverSeaOrLand)
%==========================================================================================
%--looking from the first thermal
%--breakung all the time into 10 sec pieces and looking at them 
%------we will start first piece from every thermal and not thermal and then cut off if it not full 10 sec at the end
%------the first piece will be first thermal, every thing before to the garbage
%- we will always start 1 sec after and end 2 sec before becaise it's running mean over 4 sec

%--the piece can be:
%----(1) thermal circling
%---------(a) time since started circling
%----(2) not thermal and then we record:
%---------(a) time since thermal end

%--for each piece we also get the mean height obove ground
%% ========== Loop over the thermal table and create the indexes to work with ===========================================================
AnalyseIndexes=[];
for i=1:height(ThermalTable)
 %---create indexes for pieces in thermal
    ThermalIntStartIndex=ThermalTable.IndexStart(i)+1;
    ThermalIntEndIndex=ThermalTable.IndexEnd(i)-2;
    ThermalIndexes=[[ThermalIntStartIndex:10:ThermalIntEndIndex-9];[ThermalIntStartIndex+9:10:ThermalIntEndIndex]];
    ThermalIndexes=ThermalIndexes';
    TimeSinceStartedCircling=mean(ThermalIndexes,2)-ThermalTable.IndexStart(i);
 
 %---Create the final undexes to work with
 %------[2] No thermals, flapping flight
 %------second and third coulumns are the indexes
 %------forth coulumn is the [Time Since start of the random section]
latStart=Lat(ThermalIndexes(1,1));
lonStart=Lon(ThermalIndexes(1,1));
ClimbRate=0;
Oversea=OverSeaOrLand(ThermalIndexes(1,1));
exitElevetaion=ElevationEboveG(ThermalIndexes(end,2));
TimeGlide=300; % fixed length of random sections

LT=length(ThermalIndexes(:,1));
    AnalyseIndexes=[AnalyseIndexes; [ones(LT,1)*2,ThermalIndexes,TimeSinceStartedCircling,...
        ones(LT,1)*latStart,ones(LT,1)*lonStart,ones(LT,1)*ClimbRate,ones(LT,1)*Oversea,...
        ones(LT,1)*exitElevetaion,ones(LT,1)*TimeGlide,ones(LT,1)*i]];
end

%% ========== Loop over all the indexes and calculate the mean Flap rate and Elevation ===========================================================

for f=1:length(AnalyseIndexes(:,1))
 mean_FR(f,1)=mean(FR(AnalyseIndexes(f,2):AnalyseIndexes(f,3))); 
 if mean(FR(AnalyseIndexes(f,2):AnalyseIndexes(f,3)))==0
     mean(FR(AnalyseIndexes(f,2):AnalyseIndexes(f,3)))
 end
 mean_ElevationEboveG(f,1)=nanmean(ElevationEboveG(AnalyseIndexes(f,2):AnalyseIndexes(f,3))); 
end
Answer=table(AnalyseIndexes(:,1),AnalyseIndexes(:,4),AnalyseIndexes(:,5),AnalyseIndexes(:,6),...
    AnalyseIndexes(:,7),AnalyseIndexes(:,8),AnalyseIndexes(:,9),AnalyseIndexes(:,10),mean_FR,mean_ElevationEboveG,AnalyseIndexes(:,11),...
    'VariableNames',{'CircleORNot','TimeSince','lat_start','lon_start','climb_rate_m_sec','over_sea','exit_eleG_m','time_calc_glide_sec',...
    'mean_FR','mean_ElevationEboveG','runing_ID'});
end
