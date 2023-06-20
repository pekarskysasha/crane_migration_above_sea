clear all
add=0;
load DataFor1Hz_updated1
load DataFor1Hz_updated2
% if we need to update only
add=1;
indDone=[];
if add==1
    StatsDownloaded_old=StatsDownloaded;
    load DataFor1Hz_updated3
    for sd=1:height(StatsDownloaded_old)
        indSt=find(StatsDownloaded.Individual==StatsDownloaded_old.Individual(sd) &...
            StatsDownloaded.TimestartV==StatsDownloaded_old.TimestartV(sd));
        if ~isempty(indSt) & StatsDownloaded_old.FoundData(indSt)==1
            indDone=[indDone; indSt];
        end
    end            
    DataAll=[DataAll1;DataAll2;DataAll3];
else
    DataAll=[DataAll1;DataAll2];
end 
DataAll=[DataAll1;DataAll2;DataAll3];
%% Get ACC TransConst
[nT,tT]=xlsread('../../GPS Tags\ACC_calibration\TagsParametersFinal.xlsx');
%% Load sea polygons data to know if flying over sea or over land
[Hsea, Asea] = shaperead('../../LandCover\Seas\World_Seas_IHO_v3','UseGeoCoords', true);
BS=table(Hsea(90).Lon', Hsea(90).Lat','VariableNames',{'Lon','Lat'}); %Black sea
MS_EB=table(Hsea(44).Lon', Hsea(44).Lat','VariableNames',{'Lon','Lat'}); % Mediterranean Sea - Eastern Basin
%% loop over all 1Hz
SratsT=StatsDownloaded;
SratsT(indDone,:)=[];
SratsT=SratsT(SratsT.FoundData==1,:);
Tableall=[];
TableGPSAll=[];
SensorsAll=[];
ThermalTableAll=[];
AllInfo=[];
conter=1;
StatsTCleaned=[];
for h=1:height(SratsT) % loop for all 1Hz peeces for this tag

   Ind1Hz=DataAll.HZIndex==SratsT.Ind(h);
  
    
    Data=DataAll(Ind1Hz,:);
    if isempty(Data)
        continue;
    end
    % index of tag in the TagConst
    TagConstIND=find(nT(:,1)==Data.device_id(1));
    if isempty(TagConstIND)
        TagConstIND=64; % we are using the median line + the regular orienattion of the new tags
    end
    
    %% Sensor Burst %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % exlude GPS data
    indACC=isnan(Data.Latitude);
    ACCExist=0; % variable which tels us id we have acc data
    % check if we have sensor data
    if sum(indACC)>10
        ACCExist=1;
        % Get tag constatnts
        X =(Data.acc_x(indACC)-nT(TagConstIND,2)).*nT(TagConstIND,5)*9.81*-1*(nT(TagConstIND,8)); % Calculate X-acceleration AccData m/sec2
        Y =(Data.acc_y(indACC)-nT(TagConstIND,3)).*nT(TagConstIND,6)*9.81*-1*(nT(TagConstIND,9)); % Calculate Y-acceleration AccData m/sec2
        Z =(Data.acc_z(indACC)-nT(TagConstIND,4)).*nT(TagConstIND,7)*9.81*-1*(nT(TagConstIND,10)); % Calculate Z-acceleration AccData m/sec2
        ACC_XYZ=[X,Y,Z];
        
        maldateTemp=datenum(Data.UTC_datetime);
        localmldate=maldateTemp(indACC);
        %find ACC observations
        ff=[NaN; floor(diff(localmldate)*86400)];
        indGaps=[1; find(ff>1); length(localmldate)];
        % create a time string for ACC recording begining and end
        BeginEnd_RecACC=[localmldate(indGaps(1:end-1)),localmldate(indGaps(2:end)-1)];
        Indexes=[indGaps(1:end-1),indGaps(2:end)-1];
        % Eraze short parts of busts
        IndEraze=find(abs((BeginEnd_RecACC(:,1)-BeginEnd_RecACC(:,2))*86400)<1 | ...
            abs((Indexes(:,1)-Indexes(:,2)))<15);
        BeginEnd_RecACC(IndEraze,:)=[];
        Indexes(IndEraze,:)=[];
        % create a vector for runing sensor burst number
        Counting=zeros(length(ACC_XYZ),1);
        %  convert mldate to datetime for simplicity, we will also need it laterfor the table
        VecMldate=datetime(localmldate,'ConvertFrom','datenum','format','dd-MMM-y HH:mm:ss');
        
        CalculatedVec=[];
        ACC_XYZ1=[];
        Counting=[];
        mag=[];
        VecMldate1=[];
        magTemp=[Data.mag_x(indACC),Data.mag_y(indACC),Data.mag_z(indACC)];
        if ~isempty(Indexes)  
            for rc=1:length(Indexes(:,1))
                % create time string for the sensor data based on number of samples per second
                TimeVec2=VecMldate(Indexes(rc,1):Indexes(rc,2));
                %Take the relevant data
                if rc==length(Indexes(:,1))% if the last point
                    TimeVec2=VecMldate(Indexes(rc,1):Indexes(rc,2)+1);
                    ACC_XYZ1=[ACC_XYZ1; ACC_XYZ(Indexes(rc,1):Indexes(rc,2)+1,:)];
                    mag=[mag; magTemp(Indexes(rc,1):Indexes(rc,2)+1,:)];
                    VecMldate1=[VecMldate1; VecMldate(Indexes(rc,1):Indexes(rc,2)+1)];
                else
                    ACC_XYZ1=[ACC_XYZ1; ACC_XYZ(Indexes(rc,1):Indexes(rc,2),:)];
                    mag=[mag; magTemp(Indexes(rc,1):Indexes(rc,2),:)];
                    VecMldate1=[VecMldate1; VecMldate(Indexes(rc,1):Indexes(rc,2))];
                end
                Counting=[Counting; ones(length(TimeVec2),1)*rc];
                
                Sec=[unique(TimeVec2)];
                Sec=[Sec; Sec(end)+seconds(1)];
                TimeVecDT2=[];
                for c=1:length(Sec)-1
                    numberPerSec=sum(TimeVec2==Sec(c));
                    if c==1
                        Ttemp=linspace(Sec(c)+seconds(0.5),Sec(c+1)-seconds(0.1),numberPerSec);
                        Ttemp=datetime(Ttemp','format','dd-MMM-y HH:mm:ss.SS');
                        TimeVecDT2=[TimeVecDT2; Ttemp];
                    elseif c==length(Sec)-1
                        Ttemp=linspace(Sec(c),Sec(c+1)-seconds(0.5),numberPerSec);
                        Ttemp=datetime(Ttemp','format','dd-MMM-y HH:mm:ss.SS');
                        TimeVecDT2=[TimeVecDT2; Ttemp];
                    else
                        Ttemp=linspace(Sec(c),Sec(c+1)-seconds(0.1),numberPerSec);
                        Ttemp=datetime(Ttemp','format','dd-MMM-y HH:mm:ss.SS');
                        TimeVecDT2=[TimeVecDT2; Ttemp];
                    end
                end
                CalculatedVec=[CalculatedVec; TimeVecDT2];
            end
            %% Assemble all sensors together
            VecMldate=VecMldate1;
            Tag=ones(length(VecMldate),1)*Data.device_id(1);
            HZIndex=ones(length(VecMldate),1)*SratsT.Ind(h);
            Sensors=table(Tag,HZIndex,VecMldate,CalculatedVec,ACC_XYZ1(:,1),ACC_XYZ1(:,2),ACC_XYZ1(:,3),...
                mag(:,1),mag(:,2),mag(:,3),Counting,...
                'VariableNames',{'Tag','HZIndex','TimeRaw','Time_Calculated','ACC_x','ACC_y','ACC_z','mag_x','mag_y','mag_z','Counting'});
            %% Save the Sensor data
            SensorsAll=[SensorsAll; Sensors];
        end
    end
    %% GPS (+1Hz Sensors) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    GPSdata=Data(strcmp(Data.datatype,'GPS'),:);
    %% Organize the data
    Interv=round(SratsT.MeadianInterval(h)); % What is the interval of the data
    % convert the the to sec from first point
    qt=round((datenum(GPSdata.UTC_datetime)-datenum(GPSdata.UTC_datetime(1)))*24*3600);
    % create a vector with "true" where we have data (The vector is always 1 Hz)
    qis_exist=false(max(qt)+1,1);
    qis_exist(qt+1)=true;
    % make the time string contineius and put nans in times where there are gaps
%     qtall=(0:max(qt))./(24*3600)+datenum(GPSdata.UTC_datetime(1));
%     TimeCont=datetime(qtall','ConvertFrom','datenum','format','dd-MMM-y HH:mm:ss');
    TimeCont=linspace(GPSdata.UTC_datetime(1),GPSdata.UTC_datetime(end),length(qis_exist));
    TimeCont=TimeCont';
    
    Latitude=nan(size(TimeCont));
    Longitude=nan(size(TimeCont));
    Altitude_m=nan(size(TimeCont));
    Speed_km_h=nan(size(TimeCont));
    direction_deg=nan(size(TimeCont));
    mag_x=nan(size(TimeCont));
    mag_y=nan(size(TimeCont));
    mag_z=nan(size(TimeCont));
    acc_x=nan(size(TimeCont));
    acc_y=nan(size(TimeCont));
    acc_z=nan(size(TimeCont));
    IndexCount=nan(size(TimeCont));
    %% ACC For GPS data
    XGPS =(GPSdata.acc_x-nT(TagConstIND,2)).*nT(TagConstIND,5)*9.81*-1*(nT(TagConstIND,8)); % Calculate X-acceleration AccData m/sec2
    YGPS =(GPSdata.acc_y-nT(TagConstIND,3)).*nT(TagConstIND,6)*9.81*-1*(nT(TagConstIND,9)); % Calculate Y-acceleration AccData m/sec2
    ZGPS =(GPSdata.acc_z-nT(TagConstIND,4)).*nT(TagConstIND,7)*9.81*-1*(nT(TagConstIND,10)); % Calculate Z-acceleration AccData m/sec2
% go over the vector and fill in teh data where exists
IND=1;
    for i=1:length(TimeCont) 
        if qis_exist(i) % if there is data foe this point
            Latitude(i)=GPSdata.Latitude(IND);
            Longitude(i)=GPSdata.Longitude(IND);
            Altitude_m(i)=GPSdata.Altitude_m(IND);
            Speed_km_h(i)=GPSdata.speed_km_h(IND);
            direction_deg(i)=GPSdata.direction_deg(IND);
            mag_x(i)=GPSdata.mag_x(IND);
            mag_y(i)=GPSdata.mag_y(IND);
            mag_z(i)=GPSdata.mag_y(IND);
            acc_x(i)=XGPS(IND);
            acc_y(i)=YGPS(IND);
            acc_z(i)=ZGPS(IND);
            IND=IND+1;
        else
            if ACCExist==1 % if we have ACC burst data
                IntTemp=(BeginEnd_RecACC(:,1)-datenum(TimeCont(i)))*86400;
                SI=find(IntTemp>=-10 & IntTemp<=10);
            else
                SI=[];
            end
            if ~isempty(SI)
                INDdata=find(IntTemp>=-10 & IntTemp<=10);
                if length(INDdata)>1 % if found two points, take the one after the time
                    INDdata=find(abs(IntTemp)==max(abs(IntTemp(INDdata))),1,'first');
                end
                IndexCount(i)=INDdata;
                % If we have data, fill the gaps in the Magnetometer and ACC data with the burst data
                IntervTemp=abs((datenum(Sensors.Time_Calculated)-datenum(TimeCont(i)))*86400);
                IndBurst=find(IntervTemp==min(IntervTemp));
                if length(IndBurst)>1
                    IndBurst=IndBurst(1);
                end
                % fill the data from the burst to the 1Hz table
                mag_x(i)=Sensors.mag_x(IndBurst);
                mag_y(i)=Sensors.mag_y(IndBurst);
                mag_z(i)=Sensors.mag_z(IndBurst);
                acc_x(i)=Sensors.ACC_x(IndBurst);
                acc_y(i)=Sensors.ACC_y(IndBurst);
                acc_z(i)=Sensors.ACC_z(IndBurst);
            else
                IndexCount(i)=nan;
            end
        end
    end   

    %% Check if over land or over sea
    % Black sea => 1
    % Mediterranean Sea => 2
    % Land => 0
    OverSeaOrLand=zeros(length(Longitude),1);
    IndBS=inpolygon(Longitude,Latitude,BS.Lon,BS.Lat);
    IndMS_EB=inpolygon(Longitude,Latitude,MS_EB.Lon,MS_EB.Lat);
    OverSeaOrLand(IndBS)=1;
    OverSeaOrLand(IndMS_EB)=2;
    
    %% create a table
    Tag=ones(length(TimeCont),1)*Data.device_id(1);
    HZIndex=ones(length(TimeCont),1)*SratsT.Ind(h);
    OverSeaOrLand(isnan(Latitude))=nan;
    TableGPS=table(Tag,HZIndex,TimeCont,Latitude,Longitude,Altitude_m,Speed_km_h,direction_deg,OverSeaOrLand,...
        mag_x,mag_y,mag_z,acc_x,acc_y,acc_z,IndexCount); 
%% look for missing data and exclude it
MissVec=isnan(Longitude);
MissLook=diff(MissVec);
startBurst=find(MissLook==1)+1;
endBurst=find(MissLook==-1);
IndexGaps=[startBurst,endBurst];
IndexGaps(:,3)=IndexGaps(:,2)+1-IndexGaps(:,1);
indexContinious=[[1; IndexGaps(find(IndexGaps(:,3)>10),2)+1],[IndexGaps(find(IndexGaps(:,3)>10),1);length(Longitude)]];
indexContinious(:,3)=minutes(TimeCont(indexContinious(:,2))-TimeCont(indexContinious(:,1)));
% find continious parts longer then 10 minutes
longParts=indexContinious(indexContinious(:,3)>10,:);
TableGPSfinal=[];
for lo=1:length(longParts(:,1))
    TableGPSfinal=[TableGPSfinal;  [TableGPS(longParts(lo,1):longParts(lo,2),:),...
        table(conter*ones(longParts(lo,2)-longParts(lo,1)+1,1), 'VariableNames',{'UniqueSectionCounter'})]];
    StatsTCleaned=[StatsTCleaned; [SratsT(h,14),SratsT(h,1:4),...
        table(conter,TableGPS.TimeCont(longParts(lo,1)),TableGPS.TimeCont(longParts(lo,2)),longParts(lo,3),...
        'VariableNames',{'UniqueSectionCounter','TimeStartV','TimeEndV','TimeOfBurstT'})]];
    conter=conter+1;
end
    TableGPSAll=[TableGPSAll; TableGPSfinal];
    
    
end
if add==1
TableGPSAll_2add=TableGPSAll;
SensorsAll_2add=SensorsAll;
StatsTCleaned_2add=StatsTCleaned;
load 1HzGPS_and_ACC
UnCounterAdd=max(StatsTCleaned.UniqueSectionCounter);
StatsTCleaned_2add.UniqueSectionCounter=StatsTCleaned_2add.UniqueSectionCounter+UnCounterAdd;
TableGPSAll_2add.UniqueSectionCounter=TableGPSAll_2add.UniqueSectionCounter+UnCounterAdd;
save('1HzGPS_and_ACC_2update','TableGPSAll_2add','SensorsAll_2add','StatsTCleaned_2add')
else
save('1HzGPS_and_ACC','TableGPSAll','SensorsAll','StatsTCleaned')
end
