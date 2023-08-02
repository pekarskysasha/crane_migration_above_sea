clear all
% load thermal analysis data
load Soar_new
% load bout level data
load MigrationAnalysis_Stage4_BoutLevel_L
%========================================================================
%-(A)- find sea crossing with 1HZ and check 
%---(1)- were threr thermals
%---(2)- when (a) startd migartion, (b) enterded sea and (c) crossed alt 36.13
%---(3)- how many days stayed at the stopovers (how many days are relevant to go back in anotation)
%========================================================================
%% Mediterranean Sea %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ThemalMedSea=CiclingThemalOnly(CiclingThemalOnly.OverSeaOrLand==2,:);
[a,b,c]=unique([ThemalMedSea.Animal_ID,floor(datenum(ThemalMedSea.TimeStart))],'rows');
Dates=datetime(a(:,2),'ConvertFrom','datenum','format','y-MM-dd');
Indev=a(:,1);
ThermalDaysTags=table(Indev,Dates);
[cnt_unique, unique_a] = hist(a(:,2),unique(a(:,2)));
Tag_number=cnt_unique';
Dates=datetime(unique_a,'ConvertFrom','datenum','format','y-MM-dd');
ThermalDaysTagsCount=table(Dates,Tag_number);

HZ1MedSea=TheramlStatsAll(TheramlStatsAll.OverSea==2,:);
UniqueCrane=unique(HZ1MedSea.Individual);
LostDates=[];
TimesOfMigrationFall=[];
latinsea=35.6; % point of reference in sea
for i=1:length(UniqueCrane)
    UniqueDates=unique(floor(datenum(HZ1MedSea.date(HZ1MedSea.Individual==UniqueCrane(i)))));
    for ii=1:length(UniqueDates)
        %% (1) === extract migration info ================================
        %--(1.1) - general info
        indev=UniqueCrane(i);
        Date=datetime(UniqueDates(ii),'ConvertFrom','datenum','format','y-MM-dd');
        if month(UniqueDates(ii))>9
            fall=1;
        else
            fall=0;
        end
        %--(1.2) were thermals this day? Use 1Hz data analysis
        thermal_section_num=sum(floor(datenum(HZ1MedSea.date))==UniqueDates(ii) & HZ1MedSea.Individual==UniqueCrane(i));
        mean_thermal_prop=mean(HZ1MedSea.percent_time_in_thermals(floor(datenum(HZ1MedSea.date))==UniqueDates(ii)...
            & HZ1MedSea.Individual==UniqueCrane(i)));
        if mean_thermal_prop>0
        thermals=1;
        else
         thermals=0;  
        end
        
        %% (2) === find migration bout analysis ================================
        % check for migration bout for this date & individual
        MigBout=BoutLevel_L(BoutLevel_L.animal_ID==UniqueCrane(i) &...
            floor(datenum(BoutLevel_L.bout_starting_DateTime))==UniqueDates(ii),:);
        %% (3) === check how many days stayed at the last roost ================================
        %--(3.1) extract night hours (22:23)
        RoostData=ResampledData(ResampledData.animal_ID==UniqueCrane(i) &...
            floor(datenum(ResampledData.date))>=UniqueDates(ii)-4 & floor(datenum(ResampledData.date))<UniqueDates(ii) &...
            hour(ResampledData.date_time)>=22 & hour(ResampledData.date_time)<=23,:);
        if isempty(RoostData) % No information when arrived to stopover
            Reason={'No stopover info'};
            LostDates=[LostDates; table(indev,Date,Reason)];
            continue;
        end
        % if bout found extract the first point of migration 
        if ~isempty(MigBout)
            % do we have more than one bout this day?
            if height(MigBout)>1
                Reason={'MigBout length'};
                LostDates=[LostDates; table(indev,Date,Reason)];
                continue;
            end
       %--(3.2) extract the hourly resampled data for this bout ----------------------
            %--- * it starts from the second point, add the first point
            indTemp=find(ResampledData.animal_ID==UniqueCrane(i) &...
                ResampledData.bout_num==MigBout.bout_number);            
            BoutHourlyData=ResampledData([indTemp(1); indTemp],:);
            if BoutHourlyData.lat(1)<latinsea
               Reason={'MigStart not Turkey'};
               LostDates=[LostDates; table(indev,Date,Reason)];
               continue;
            end
        else 
            Reason={'Missing MigBout'};
            LostDates=[LostDates; table(indev,Date,Reason)];           
%             DateHourlyData=ResampledData(ResampledData.animal_ID==UniqueCrane(i) &...
%                 floor(datenum(ResampledData.date))==UniqueDates(ii),:);
            continue;
        end
        %--(3.3) calculate distance between roost and first migration bout point for 5 days preior bout ----------------------
        DatesRoost=unique(RoostData.date);
        Distance2travelPoint=zeros(length(DatesRoost),1);
        TableDistnaces=table(DatesRoost,Distance2travelPoint);
        for d=1:length(DatesRoost)
            %--(3.3.1) calculte the roost location
            roost_lon=mean(RoostData.lon(RoostData.date==DatesRoost(d)));
            roost_lat=mean(RoostData.lat(RoostData.date==DatesRoost(d)));
            %--(3.3.2) find the distance from the roosing points to the first migration bout point
            TableDistnaces.Distance2travelPoint(d)=distance(roost_lat,roost_lon,...
                BoutHourlyData.lat(1),BoutHourlyData.lon(1),wgs84Ellipsoid)/1000;
            
        end
        %--(3.3.3) find the last day that the roost was more then 100 km away from the roost of the last day
            First_day_this_roost=find(TableDistnaces.Distance2travelPoint>TableDistnaces.Distance2travelPoint(end)+100,1,'last');
            if ~isempty(First_day_this_roost)
                Days_at_stopover=height(TableDistnaces)-First_day_this_roost-1;
            else
                % put maximal days checked
                Days_at_stopover=3;
            end
        %--(3.3.4) fing the heading of the bird between the firt point and the nearest point to latitide of insea
        inPoint4ang=find(abs(BoutHourlyData.lat-latinsea)==min(abs(BoutHourlyData.lat-latinsea)),1,'first');
        ang=atan2(BoutHourlyData.lat(inPoint4ang)-BoutHourlyData.lat(1),BoutHourlyData.lon(inPoint4ang)-BoutHourlyData.lon(1));
        %% (4) === find the relevant bout points for anotation ================================ 
        % find the first point before migration bout began
        PointBeforeLeft=ResampledData(find(ResampledData.animal_ID==UniqueCrane(i) &...
            ResampledData.date_time==MigBout.bout_starting_DateTime)-1,:);
        datetime_start_migration_utc=PointBeforeLeft.date_time(1);
        datetime_enter_sea_utc=BoutHourlyData.date_time(find(BoutHourlyData.sea_crossing==2,1,'first')+1);
        datetime_cross_3613_lat_utc=BoutHourlyData.date_time(find(BoutHourlyData.lat<36.13,1,'first')+1);
        %% (5) === fill table to make anotation later ================================ 
        lat_start_migration=PointBeforeLeft.lat(1);
        lon_start_migration=PointBeforeLeft.lon(1);
        TimesOfMigrationFall=[TimesOfMigrationFall; table(indev,Date,fall,thermal_section_num,mean_thermal_prop,thermals,...
                                  lat_start_migration,lon_start_migration,datetime_start_migration_utc,...
                                  datetime_enter_sea_utc,datetime_cross_3613_lat_utc,Days_at_stopover,ang)];
    end
end
TimesOfMigrationFall=TimesOfMigrationFall(TimesOfMigrationFall.fall==1,:);
[a,b,c]=unique(TimesOfMigrationFall(:,[2,6]),'rows');
DatesMig=unique(a.Date);
MigDaysTheraml=[];
for j=1:length(DatesMig)
    Th=a.thermals(a.Date==DatesMig(j));
    if sum(Th==1)>0
        Thermals=1;
        Date_mig=DatesMig(j);
        MigDaysTheraml=[MigDaysTheraml; table(Date_mig,Thermals)];
    else
        Thermals=0;
        Date_mig=DatesMig(j);
        MigDaysTheraml=[MigDaysTheraml; table(Date_mig,Thermals)];
    end 
end
writetable(MigDaysTheraml,'MigDaysTheraml.csv')
%% --- create table for track annotation---------------------
departure=[];
seaEntrance=[];
Insea=[];
for i=1:height(TimesOfMigrationFall)
    % how many days stayed in this stopover
    daysBackLoop=[-TimesOfMigrationFall.Days_at_stopover(1):1:1]';
    daysBackLoop=sortrows(daysBackLoop,'descend');
    for ii=1:length(daysBackLoop)
        daysBack = daysBackLoop(ii);
        thermap_presence=TimesOfMigrationFall.thermals(i);
        %-- departure
        dateTemp=datenum(TimesOfMigrationFall.datetime_start_migration_utc(i))+daysBackLoop(ii);
        datetime_utc=datetime(dateTemp,'ConvertFrom','datenum','format','y-MM-dd HH:mm:ss');
        lon=TimesOfMigrationFall.lon_start_migration(i);
        lat=TimesOfMigrationFall.lat_start_migration(i);
        mig_ang=TimesOfMigrationFall.ang(i);
        departure=[departure; [TimesOfMigrationFall(i,1:3),table(thermap_presence,daysBack,lat,lon,mig_ang,datetime_utc)]];
        
        %-- seaEntrance
        dateTemp=datenum(TimesOfMigrationFall.datetime_enter_sea_utc(i))+daysBackLoop(ii);
        datetime_utc=datetime(dateTemp,'ConvertFrom','datenum','format','y-MM-dd HH:mm:ss');
        if TimesOfMigrationFall.lon_start_migration(i)>34.5 % adana
            lon=35.5;
            lat=36.4;
        elseif TimesOfMigrationFall.lat_start_migration(i)>36.5 % tuz
            lon=34.1;
            lat=36.15;
        else % cyprus
            lon=33.8;
            lat=34.75;
        end
        if ~isnan(datenum(TimesOfMigrationFall.datetime_enter_sea_utc(i)))
            seaEntrance=[seaEntrance; [TimesOfMigrationFall(i,1:3),table(thermap_presence,daysBack,lat,lon,mig_ang,datetime_utc)]];
        end
        
        %-- In sea
        dateTemp=datenum(TimesOfMigrationFall.datetime_cross_3613_lat_utc(i))+daysBackLoop(ii);
        datetime_utc=datetime(dateTemp,'ConvertFrom','datenum','format','y-MM-dd HH:mm:ss');
        lon=35.3;
        lat=35.6;
        if ~isnan(datenum(TimesOfMigrationFall.datetime_cross_3613_lat_utc(i)))
        Insea=[Insea; [TimesOfMigrationFall(i,1:3),table(thermap_presence,daysBack,lat,lon,mig_ang,datetime_utc)]];
        end
    end
end

Departure_MedSea=departure;
SeaEntrance_MedSea=seaEntrance;
Insea_MedSea=Insea;

figure
plot(departure.lon,departure.lat,'.','Color','r','MarkerSize',12)
hold on
plot(seaEntrance.lon,seaEntrance.lat,'.','Color','y','MarkerSize',12)
plot(Insea.lon,Insea.lat,'.','Color','g','MarkerSize',12)
plot_google_map('MapType','satellite')
legend('departure','sea entrance', 'in sea')
%% Black Sea %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ThemalBKSea=CiclingThemalOnly(CiclingThemalOnly.OverSeaOrLand==1,:);
[a,b,c]=unique([ThemalBKSea.Animal_ID,floor(datenum(ThemalBKSea.TimeStart))],'rows');
Dates=datetime(a(:,2),'ConvertFrom','datenum','format','y-MM-dd');
Indev=a(:,1);
ThermalDaysTags=table(Indev,Dates);

HZ1BKSea=TheramlStatsAll(TheramlStatsAll.OverSea==1,:);
UniqueCrane=unique(HZ1BKSea.Individual);
LostDates=[];
TimesOfMigrationFall=[];
for i=1:length(UniqueCrane)
    UniqueDates=unique(floor(datenum(HZ1BKSea.date(HZ1BKSea.Individual==UniqueCrane(i)))));
    for ii=1:length(UniqueDates)
        %% (1) === extract migration info ================================
        %--(1.1) - general info
        indev=UniqueCrane(i);
        Date=datetime(UniqueDates(ii),'ConvertFrom','datenum','format','y-MM-dd');
        if month(UniqueDates(ii))>9
            fall=1;
        else
            fall=0;
        end
        %--(1.2) were thermals this day? Use 1Hz data analysis
        thermal_section_num=sum(floor(datenum(HZ1BKSea.date))==UniqueDates(ii) & HZ1BKSea.Individual==UniqueCrane(i));
        mean_thermal_prop=mean(HZ1BKSea.percent_time_in_thermals(floor(datenum(HZ1BKSea.date))==UniqueDates(ii)...
            & HZ1BKSea.Individual==UniqueCrane(i)));
        if mean_thermal_prop>0
        thermals=1;
        else
         thermals=0;  
        end
        
        %% (2) === find migration bout bout analysis ================================
        % check for migration bout for this date & individual
        MigBout=BoutLevel_L(BoutLevel_L.animal_ID==UniqueCrane(i) &...
            floor(datenum(BoutLevel_L.bout_starting_DateTime))==UniqueDates(ii),:);
        %% (3) === check how many days stayed at the last roost ================================
        %--(3.1) extract night hours (22:23)
        RoostData=ResampledData(ResampledData.animal_ID==UniqueCrane(i) &...
            floor(datenum(ResampledData.date))>=UniqueDates(ii)-4 & floor(datenum(ResampledData.date))<UniqueDates(ii) &...
            hour(ResampledData.date_time)>=22 & hour(ResampledData.date_time)<=23,:);
        if isempty(RoostData) % No information when arrived to stopover
            Reason={'No stopover info'};
            LostDates=[LostDates; table(indev,Date,Reason)];
            continue;
        end
        % if bout found extract the first point of migration 
        if ~isempty(MigBout)
            % do we have more than one bout this day?
            if height(MigBout)>1
                Reason={'MigBout length'};
                LostDates=[LostDates; table(indev,Date,Reason)];
                continue;
            end
       %--(3.2) extract the hourly resamples data for this bout ----------------------
            %--- * it starts from the second point, add the first point
            indTemp=find(ResampledData.animal_ID==UniqueCrane(i) &...
                ResampledData.bout_num==MigBout.bout_number);            
            BoutHourlyData=ResampledData([indTemp(1); indTemp],:);
            if sum(BoutHourlyData.sea_crossing==1)==0
                Reason={'No actual migration over sea'};
                LostDates=[LostDates; table(indev,Date,Reason)];
                continue;
            end
        else 
            Reason={'Missing MigBout'};
            LostDates=[LostDates; table(indev,Date,Reason)];           
%             DateHourlyData=ResampledData(ResampledData.animal_ID==UniqueCrane(i) &...
%                 floor(datenum(ResampledData.date))==UniqueDates(ii),:);
            continue;
        end
        %--(3.3) calculate distance between roost and first migration bout point for 5 days preior bout ----------------------
        DatesRoost=unique(RoostData.date);
        Distance2travelPoint=zeros(length(DatesRoost),1);
        TableDistnaces=table(DatesRoost,Distance2travelPoint);
        for d=1:length(DatesRoost)
            %--(3.3.1) calculte the roost location
            roost_lon=mean(RoostData.lon(RoostData.date==DatesRoost(d)));
            roost_lat=mean(RoostData.lat(RoostData.date==DatesRoost(d)));
            %--(3.3.2) find the distance from the roosing points to the first migration bout point
            TableDistnaces.Distance2travelPoint(d)=distance(roost_lat,roost_lon,...
                BoutHourlyData.lat(1),BoutHourlyData.lon(1),wgs84Ellipsoid)/1000;
            
        end
        %--(3.3.3) find the last day that the roost was more then 100 km away from the roost of the last day
            First_day_this_roost=find(TableDistnaces.Distance2travelPoint>TableDistnaces.Distance2travelPoint(end)+100,1,'last');
            if ~isempty(First_day_this_roost)
                Days_at_stopover=height(TableDistnaces)-First_day_this_roost-1;
            else
                % put maximal days checked
                Days_at_stopover=3;
            end
        %% (4) === find the relevant bout points for anotation ================================ 
        % find the first point before migration bout began
        PointBeforeLeft=ResampledData(find(ResampledData.animal_ID==UniqueCrane(i) &...
            ResampledData.date_time==MigBout.bout_starting_DateTime)-1,:);
        datetime_start_migration_utc=PointBeforeLeft.date_time(1);
        datetime_enter_sea_utc=BoutHourlyData.date_time(find(BoutHourlyData.sea_crossing==1,1,'first')+1);
        datetime_cross_435_lat_utc=BoutHourlyData.date_time(find(BoutHourlyData.lat<43.5,1,'first')+1);
        %% (5) === fill table to make anotation later ================================ 
        lat_start_migration=PointBeforeLeft.lat(1);
        lon_start_migration=PointBeforeLeft.lon(1);
        lat_SeaEntrance=BoutHourlyData.lat(find(BoutHourlyData.sea_crossing==1,1,'first')+1);
        lon_SeaEntrance=BoutHourlyData.lon(find(BoutHourlyData.sea_crossing==1,1,'first')+1);
        TimesOfMigrationFall=[TimesOfMigrationFall; table(indev,Date,fall,thermal_section_num,mean_thermal_prop,thermals,...
                                  lat_start_migration,lon_start_migration,datetime_start_migration_utc,...
                                  datetime_enter_sea_utc,lat_SeaEntrance,lon_SeaEntrance,...
                                  datetime_cross_435_lat_utc,Days_at_stopover)];
    end
end
TimesOfMigrationFall=TimesOfMigrationFall(TimesOfMigrationFall.fall==1,:);
%% --- create table for track annotation---------------------
departure=[];
seaEntrance=[];
Insea=[];
for i=1:height(TimesOfMigrationFall)
    % how many days stayed in this stopover
    Days_at_stopover=[-TimesOfMigrationFall.Days_at_stopover(1):1:1]';
    Days_at_stopover=sortrows(Days_at_stopover,'descend');
    for ii=1:length(daysBackLoop)
        daysBack = daysBackLoop(ii);
        thermap_presence=TimesOfMigrationFall.thermals(i);
        %-- departure
        dateTemp=datenum(TimesOfMigrationFall.datetime_start_migration_utc(i))+daysBackLoop(ii);
        datetime_utc=datetime(dateTemp,'ConvertFrom','datenum','format','y-MM-dd HH:mm:ss');
        lon=TimesOfMigrationFall.lon_start_migration(i);
        lat=TimesOfMigrationFall.lat_start_migration(i);
        departure=[departure; [TimesOfMigrationFall(i,1:3),table(thermap_presence,daysBack,lat,lon,datetime_utc)]];
        
        %-- seaEntrance
            lon= TimesOfMigrationFall.lon_SeaEntrance(i);
            lat= TimesOfMigrationFall.lat_SeaEntrance(i);
     
        if ~isnan(datenum(TimesOfMigrationFall.datetime_enter_sea_utc(i)))
            seaEntrance=[seaEntrance; [TimesOfMigrationFall(i,1:3),table(thermap_presence,daysBack,lat,lon,datetime_utc)]];
        end
        
        %-- In sea
        dateTemp=datenum(TimesOfMigrationFall.datetime_cross_435_lat_utc(i))+daysBackLoop(ii);
        datetime_utc=datetime(dateTemp,'ConvertFrom','datenum','format','y-MM-dd HH:mm:ss');
        lon=33.98;
        lat=43.5;
        if ~isnan(datenum(TimesOfMigrationFall.datetime_cross_435_lat_utc(i)))
        Insea=[Insea; [TimesOfMigrationFall(i,1:3),table(thermap_presence,daysBack,lat,lon,datetime_utc)]];
        end
    end
end

Departure_BkSea=departure;
SeaEntrance_BkSea=seaEntrance;
Insea_BkSea=Insea;


figure
plot(departure.lon,departure.lat,'.','Color','r','MarkerSize',12)
hold on
plot(seaEntrance.lon,seaEntrance.lat,'.','Color','y','MarkerSize',12)
plot(Insea.lon,Insea.lat,'.','Color','g','MarkerSize',12)
plot_google_map('MapType','satellite')
legend('departure','sea entrance', 'in sea')

save('TimeScales','Departure_MedSea','SeaEntrance_MedSea','Insea_MedSea','Departure_BkSea','SeaEntrance_BkSea','Insea_BkSea')

%% load annotated
load 'ThermalAnalysis/Working folder/TimeScales_annotated'
%- Join tables togeher
Departure_BkSea.time_point(:,1)=1;
SeaEntrance_BkSea.time_point(:,1)=2;
Insea_BkSea.time_point(:,1)=3;
AnnotatedTimePointsBKSea=[Departure_BkSea;SeaEntrance_BkSea;Insea_BkSea];

%--Med sea
Departure_MedSea.time_point(:,1)=1;
SeaEntrance_MedSea.time_point(:,1)=2;
Insea_MedSea.time_point(:,1)=3;
%=============================================================================================
%--if using the old anotaton add the angel data
Departure_MedSea(Departure_MedSea.indev==170991 & Departure_MedSea.Date=='2019-10-25',:)=[];
Departure_MedSea=[Departure_MedSea(:,1:7),Departure_MedSea1(:,8),Departure_MedSea(:,8:31)];
Insea_MedSea(Insea_MedSea.indev==170991 & Insea_MedSea.Date=='2019-10-25',:)=[];
Insea_MedSea=[Insea_MedSea(:,1:7),Insea_MedSea1(:,8),Insea_MedSea(:,8:31)];
SeaEntrance_MedSea(SeaEntrance_MedSea.indev==170991 & SeaEntrance_MedSea.Date=='2019-10-25',:)=[];
SeaEntrance_MedSea=[SeaEntrance_MedSea(:,1:7),SeaEntrance_MedSea1(:,8),SeaEntrance_MedSea(:,8:31)];
%=============================================================================================
AnnotatedTimePointsMedSea=[Departure_MedSea;SeaEntrance_MedSea;Insea_MedSea];

%-- calculate tail wind-------
% V => Velocity of the north-south (meridoinal) component of wind. Positive values indicate south to north flow (m/s) => Y
% U => Velocity of the east-west (zonal) component of wind. Positive values indicate west to east flow (m/s) => X
wang=atan2(AnnotatedTimePointsMedSea.v925,AnnotatedTimePointsMedSea.u925);
wvel=hypot(AnnotatedTimePointsMedSea.u925,AnnotatedTimePointsMedSea.v925);
winddir = atan2d(AnnotatedTimePointsMedSea.v925,AnnotatedTimePointsMedSea.u925) + 360*(AnnotatedTimePointsMedSea.v925<0);
% tail wind
tw=cos(wang-AnnotatedTimePointsMedSea.mig_ang).*wvel;
% side wind
sw=abs(sin(wang-ang)).*wvel;
% add wind 
AnnotatedTimePointsMedSea.tw=tw;
AnnotatedTimePointsMedSea.sw=sw;

writetable(AnnotatedTimePointsBKSea,'AnnotatedTimePointsBKSea.csv')
writetable(AnnotatedTimePointsMedSea,'AnnotatedTimePointsMedSea.csv')