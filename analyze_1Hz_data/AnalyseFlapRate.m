function [Answer]=AnalyseFlapRate(ThermalTable,Lat,Lon,T,ElevationEboveG,FR,OverSeaOrLand)
%==========================================================================================
%--looking from the first thermal
%--breaking all the time into 10 sec pieces and looking at them 
%------we will start first piece from every thermal and not thermal and then cut off if it not full 10 sec at the end
%------the first piece will be first thermal, every thing before to the garbage
%- we will always start 1 sec after and end 2 sec before because the flap rate is running mean over 4 sec

%--the section (length 10 src.) can be:
%----(1) thermal circling
%---------(a) time since started circling
%----(2) not thermal and then we record:
%---------(a) time since circling ended end

%--for each piece we also get the mean height obove ground
%% ========== Loop over the thermal table and create the indexes to work with ===========================================================
Answer=[];
AnalyseIndexes=[];
for i=1:height(ThermalTable)
 %---create indexes for pieces in thermal
    ThermalIntStartIndex=ThermalTable.IndexStart(i)+1;
    ThermalIntEndIndex=ThermalTable.IndexEnd(i)-2;
    ThermalIndexes=[[ThermalIntStartIndex:10:ThermalIntEndIndex-9];[ThermalIntStartIndex+9:10:ThermalIntEndIndex]];
    ThermalIndexes=ThermalIndexes';
    TimeSinceStartedCircling=mean(ThermalIndexes,2)-ThermalTable.IndexStart(i);
 %---create indexes for pieces not in thermal 
 %--------between the end of the thermal and begining of next themal
    if i<height(ThermalTable)
        NotThermalIntStartIndex=ThermalTable.IndexEnd(i)+1;
        NotThermalIntEndIndex=ThermalTable.IndexStart(i+1)-3;
    else %if it the last thermal look until the end
        NotThermalIntStartIndex=ThermalTable.IndexEnd(i)+2;
        NotThermalIntEndIndex=length(Lat)-2;
    end
    NotThermalIndexes=[[NotThermalIntStartIndex:10:NotThermalIntEndIndex-9];[NotThermalIntStartIndex+9:10:NotThermalIntEndIndex]];
    NotThermalIndexes=NotThermalIndexes';
    TimeSinceTheraml=mean(NotThermalIndexes,2)-ThermalTable.IndexEnd(i)+1;
 
 %-- Skip if:
   %- soar length < 1 min
   %- glide length < 1/2 min or > 15 min
   if isempty(NotThermalIndexes) | isempty(ThermalIndexes)
       continue
   end
 glidelength=NotThermalIndexes(end,2)-NotThermalIndexes(1,1);
 soarlength=ThermalIndexes(end,2)-ThermalIndexes(1,1);
 if soarlength<60 | glidelength<30 | glidelength>15*60
     continue
 end

 %---Create the final indexes to work with
 %------[1] is thermal
 %------[0] is not thermal
 %------second and third coulumns are the indexes
 %------forth coulumn is the [Time Since Theraml ended] for NOT thermal and [Time Since Circling Started] for thermals
latStart=Lat(ThermalIndexes(1,1));
lonStart=Lon(ThermalIndexes(1,1));
ClimbRate=ThermalTable.ClimbRate(i);
Oversea=OverSeaOrLand(ThermalIndexes(1,1));
exitElevetaion=ElevationEboveG(ThermalIndexes(end,2));
TimeGlide=NotThermalIndexes(end,2)-NotThermalIndexes(1,1);

LT=length(ThermalIndexes(:,1));
LNT=length(NotThermalIndexes(:,1));  

    AnalyseIndexes=[AnalyseIndexes; [ones(LT,1),ThermalIndexes,TimeSinceStartedCircling,...
        ones(LT,1)*latStart,ones(LT,1)*lonStart,ones(LT,1)*ClimbRate,ones(LT,1)*Oversea,ones(LT,1)*exitElevetaion,...
        ones(LT,1)*TimeGlide,ones(LT,1)*i]];
    AnalyseIndexes=[AnalyseIndexes; [zeros(LNT,1),NotThermalIndexes,TimeSinceTheraml,...
        ones(LNT,1)*latStart,ones(LNT,1)*lonStart,ones(LNT,1)*ClimbRate,ones(LNT,1)*Oversea,...
        ones(LNT,1)*exitElevetaion,ones(LNT,1)*TimeGlide,ones(LNT,1)*i]];
end

%% ========== Loop over all the indexes and calculate the mean Flap rate and Elevation ===========================================================
if isempty(AnalyseIndexes)
    return
end
for f=1:length(AnalyseIndexes(:,1))
 mean_FR(f,1)=mean(FR(AnalyseIndexes(f,2):AnalyseIndexes(f,3))); 
 mean_ElevationEboveG(f,1)=nanmean(ElevationEboveG(AnalyseIndexes(f,2):AnalyseIndexes(f,3))); 
end
Answer=table(AnalyseIndexes(:,1),AnalyseIndexes(:,4),AnalyseIndexes(:,5),AnalyseIndexes(:,6),...
    AnalyseIndexes(:,7),AnalyseIndexes(:,8),AnalyseIndexes(:,9),AnalyseIndexes(:,10),mean_FR,mean_ElevationEboveG,AnalyseIndexes(:,11),...
    'VariableNames',{'CircleORNot','TimeSince','lat_start','lon_start','climb_rate_m_sec','over_sea','exit_eleG_m','time_calc_glide_sec',...
    'mean_FR','mean_ElevationEboveG','runing_ID'});
end
