function [ThermalA]=equal_interval_thermal(Individual,ThermalTable,TableGPSCont)
%--Analyses the statistics of thermals, but only if this travel part is longer then 10 minutes
ThermalA=[];
LocationInd=[];
lo=0;
%% --- devision into 10 minutes sections -------------------
while lo<height(TableGPSCont)
    LocationInd=[LocationInd;[lo+1,lo+601]];
    lo=lo+601;
end
% if last index exceeds data length
if LocationInd(end,2)>height(TableGPSCont)
    LocationInd(end,2)=height(TableGPSCont);
end
%----Loop over all the indexes of land/sea
for i=1:length(LocationInd(:,1))
    date=datetime(TableGPSCont.TimeCont(LocationInd(i,1)),'format','dd-MMM-y HH:mm:ss');
    time_of_part_min=minutes(TableGPSCont.TimeCont(LocationInd(i,2))-TableGPSCont.TimeCont(LocationInd(i,1)));
    if time_of_part_min<10
        continue;
    end
    %-- (a) flight parameters
    vg=hypot(TableGPSCont.Nir_X(LocationInd(i,2))-TableGPSCont.Nir_X(LocationInd(i,1)),TableGPSCont.Nir_Y(LocationInd(i,2))-TableGPSCont.Nir_Y(LocationInd(i,1)))./(time_of_part_min*60);
    ang=atan2(TableGPSCont.Nir_Y(LocationInd(i,2))-TableGPSCont.Nir_Y(LocationInd(i,1)),TableGPSCont.Nir_X(LocationInd(i,2))-TableGPSCont.Nir_X(LocationInd(i,1)));
    wu=mean(TableGPSCont.u925(LocationInd(i,1):LocationInd(i,2)));
    wv=mean(TableGPSCont.v925(LocationInd(i,1):LocationInd(i,2)));
    
    % atan2(Y,X)
    % V => Velocity of the north-south (meridoinal) component of wind. Positive values indicate south to north flow (m/s) => Y
    % U => Velocity of the east-west (zonal) component of wind. Positive values indicate west to east flow (m/s) => X
    wang=atan2(wv,wu);
    wvel=hypot(wu,wv);
    tw=cos(wang-ang).*wvel;
    sw=abs(sin(wang-ang)).*wvel;
    va=hypot(vg-tw,sw);
    
    TotalDistance=sum(distance(TableGPSCont.Interpolated_Lat(LocationInd(i,1):(LocationInd(i,2)-1)), TableGPSCont.Interpolated_Lon(LocationInd(i,1):(LocationInd(i,2)-1)),...
        TableGPSCont.Interpolated_Lat(LocationInd(i,1)+1:(LocationInd(i,2))),TableGPSCont.Interpolated_Lon(LocationInd(i,1)+1:(LocationInd(i,2))),wgs84Ellipsoid))/1000; %km
    %-- (b) over sea or land
    % check if the section has both land and sea, and if it has, drop it
    NonNanData=TableGPSCont.OverSeaOrLand(LocationInd(i,1):(LocationInd(i,2)));
    NonNanData=NonNanData(~isnan(NonNanData));
    sealand=unique(NonNanData);
    if length(sealand)>1
        continue;
    end
    OverSea=nanmedian(TableGPSCont.OverSeaOrLand(LocationInd(i,1):(LocationInd(i,2))));
    %-- (b) enviromental parameters
    Mean_blh=nanmean(TableGPSCont.blh(LocationInd(i,1):LocationInd(i,2)));
    Mean_msl=nanmean(TableGPSCont.msl(LocationInd(i,1):LocationInd(i,2)));
    Mean_DeltaT=nanmean(TableGPSCont.sstDiff(LocationInd(i,1):LocationInd(i,2)));
    meanflapRate=nanmean(TableGPSCont.running_Flap_rate_4sec(LocationInd(i,1):LocationInd(i,2)));
    
    %--(c) Time of day (if any point during night, define as night)
    if sum(TableGPSCont.day_light(LocationInd(i,1):LocationInd(i,2))==0)>0
        day_time=0;
    else
        day_time=1;
    end
    %--(d) section start after sunrize
    [SunRiseSet]=suncycle(TableGPSCont.Interpolated_Lat(LocationInd(i,1)),...
        TableGPSCont.Interpolated_Lon(LocationInd(i,1)),datenum(TableGPSCont.TimeCont(LocationInd(i,1))));
    SunRiseSet2use=SunRiseSet/24; %sunrise and sunset in UTC
    sunriseThisDay=SunRiseSet2use(1)+floor(datenum(TableGPSCont.TimeCont(LocationInd(i,1))));
    time_sice_sunrize_h=(datenum(TableGPSCont.TimeCont(LocationInd(i,1)))-sunriseThisDay)*86400/60/60;
    sunsetThisDay=SunRiseSet2use(2)+floor(datenum(TableGPSCont.TimeCont(LocationInd(i,1))));
    time_sice_sunset_h=(datenum(TableGPSCont.TimeCont(LocationInd(i,1)))-sunsetThisDay)*86400/60/60;
    %----- if section after midnight but before sunrise, calculte for previous fay
    if time_sice_sunrize_h<0 & time_sice_sunset_h<0
        [SunRiseSet]=suncycle(TableGPSCont.Interpolated_Lat(LocationInd(i,1)),...
            TableGPSCont.Interpolated_Lon(LocationInd(i,1)),datenum(TableGPSCont.TimeCont(LocationInd(i,1)))-1);
        SunRiseSet2use=SunRiseSet/24; %sunrise and sunset in UTC
        sunsetThisDay=SunRiseSet2use(2)+floor(datenum(TableGPSCont.TimeCont(LocationInd(i,1))))-1;
        time_sice_sunset_h=(datenum(TableGPSCont.TimeCont(LocationInd(i,1)))-sunsetThisDay)*86400/60/60;
    end
    %--(e) thermal statistics
    if ~isempty(ThermalTable) % check if we even have a thermal table
        ThermalTablePart=ThermalTable(ThermalTable.IndexStart>=LocationInd(i,1) & ThermalTable.IndexEnd<LocationInd(i,2),:);
    else
        ThermalTablePart=[];
    end
    
    if ~isempty(ThermalTablePart)
        percent_time_in_thermals=sum(ThermalTablePart.TimeInThermal)/(time_of_part_min*60);
        time_in_thermals_sec=sum(ThermalTablePart.TimeInThermal);
        mean_thermal_length_sec=mean(ThermalTablePart.TimeInThermal);
        mean_climb_rate=mean(ThermalTablePart.ClimbRate);
        number_of_thermals=height(ThermalTablePart);
        % average Maximal elevetion above terrain in thermals
        maxEG=[];
        for ii=1:height(ThermalTablePart)
            maxEG=[maxEG; nanmax(TableGPSCont.AleAboveTerrain(ThermalTablePart.IndexStart(ii):ThermalTablePart.IndexEnd(ii)))];
        end
        mean_max_elvation_ebove_ground=mean(maxEG);
    else
        percent_time_in_thermals=0;
        time_in_thermals_sec=0;
        mean_thermal_length_sec=nan;
        mean_climb_rate=nan;
        mean_max_elvation_ebove_ground=nan;
        number_of_thermals=0;
    end
    UniqueSectionCounter=TableGPSCont.UniqueSectionCounter(1);
    lon_start=TableGPSCont.Interpolated_Lon(LocationInd(i,1));
    lat_start=TableGPSCont.Interpolated_Lat(LocationInd(i,1));
    ThermalA=[ThermalA;...
        table(Individual,UniqueSectionCounter,date,lat_start,lon_start,...
        OverSea,day_time,time_sice_sunrize_h,time_sice_sunset_h,time_of_part_min,TotalDistance,vg,va,...
        tw,sw,meanflapRate,percent_time_in_thermals,time_in_thermals_sec,...
        mean_thermal_length_sec,mean_climb_rate,mean_max_elvation_ebove_ground,...
        Mean_blh,Mean_msl,Mean_DeltaT,wvel,wang,number_of_thermals)] ;
end
end