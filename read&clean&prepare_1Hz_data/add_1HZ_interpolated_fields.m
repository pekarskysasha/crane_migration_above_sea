clear all
add=0;
% if we need to update only
add=1;
if add==0
    load 1HzGPS_and_ACC
else
    load 1HzGPS_and_ACC_2update
    TableGPSAll=TableGPSAll_2add;
    SensorsAll=SensorsAll_2add;
    StatsTCleaned=StatsTCleaned_2add;
end
%% Add some missing fields
TableGPSAll.CountIndex_SensorData=TableGPSAll.IndexCount;
TableGPSAll.Index_SensorData=TableGPSAll.HZIndex;
%% interpolate and add x,y
tracks=unique(TableGPSAll.UniqueSectionCounter);
%- loop per track
for i=1:length(tracks)
    IND1=TableGPSAll.UniqueSectionCounter==tracks(i);
    IND2=TableGPSAll.UniqueSectionCounter==tracks(i) & ~isnan(TableGPSAll.Latitude); % withouyt nanas to interpolate
    mldateInt=datenum(TableGPSAll.TimeCont(IND1));
    mldate=datenum(TableGPSAll.TimeCont(IND2));
    lat=TableGPSAll.Latitude(IND2);
    lon=TableGPSAll.Longitude(IND2);
    ele=TableGPSAll.Altitude_m(IND2);
    %% ----  interpolate ---------------------------------------
    DataInterp=interp1(mldate,[lon,lat,ele],mldateInt);
    latInt=DataInterp(:,2);
    lonInt=DataInterp(:,1);
    eleInt=DataInterp(:,3);
    % calculate x, y
    [d,a]=distance(latInt(1),lonInt(1),latInt,lonInt,wgs84Ellipsoid);
    x=sin(deg2rad(a)).*d; y=cos(deg2rad(a)).*d;
    TableGPSAll.Interpolated_Lat(IND1,1)=latInt;
    TableGPSAll.Interpolated_Lon(IND1,1)=lonInt;
    TableGPSAll.InterpolatedElevation(IND1,1)=eleInt;
    TableGPSAll.Nir_X(IND1,1)=x;
    TableGPSAll.Nir_Y(IND1,1)=y;

end
%% save
TableGPSAllNew=TableGPSAll;
if add==0
   save('1HzGPS_and_ACC_final','TableGPSAllNew','SensorsAll','StatsTCleaned')
else
    TableGPSAllNew_2add=TableGPSAllNew;
    SensorsAll_2add=SensorsAll;
    StatsTCleaned_2add=StatsTCleaned;
    save('1HzGPS_and_ACC_2update','TableGPSAllNew_2add','SensorsAll_2add','StatsTCleaned_2add')
end


%% --- elevation above ground -----------------
        % -----------For nev-Data Movebank-------------------------------------
mladateTV=datetime(datenum(TableGPSAll.TimeCont),'ConvertFrom','datenum','format','y-MM-dd HH:mm:ss.SSS');
Tabl=[TableGPSAll(:,[1,17]), ...
    table(mladateTV,TableGPSAll.Interpolated_Lat,TableGPSAll.Interpolated_Lon,TableGPSAll.InterpolatedElevation,...
    'VariableNames',{'timestamp','location_lat','location_long','height_above_msl'})];
%%- devision A
Tabl1=Tabl(1:564589,:);
Tabl(1:height(Tabl1),:)=[];
Devision=height(Tabl)/5;
Tabl2=Tabl(1:Devision,:);
Tabl(1:Devision,:)=[];
Tabl3=Tabl(1:Devision,:);
Tabl(1:Devision,:)=[];
Tabl4=Tabl(1:Devision,:);
Tabl(1:Devision,:)=[];
Tabl5=Tabl(1:Devision,:);
Tabl(1:Devision,:)=[];
Tabl6=Tabl(1:Devision,:);

writetable(Tabl1,'TableGPSAll_1.txt','Delimiter',',')
writetable(Tabl2,'TableGPSAll_2.txt','Delimiter',',')
writetable(Tabl3,'TableGPSAll_3.txt','Delimiter',',')
writetable(Tabl4,'TableGPSAll_4.txt','Delimiter',',')
writetable(Tabl5,'TableGPSAll_5.txt','Delimiter',',')
writetable(Tabl6,'TableGPSAll_6.txt','Delimiter',',')

%%- devision B
Devision=height(Tabl)/4;
Tabl1=Tabl(1:Devision,:);
Tabl(1:Devision,:)=[];
Tabl2=Tabl(1:Devision,:);
Tabl(1:Devision,:)=[];
Tabl3=Tabl(1:Devision,:);
Tabl(1:Devision,:)=[];
Tabl4=Tabl(1:Devision,:);

writetable(Tabl1,'TableGPSAll_7.txt','Delimiter',',')
writetable(Tabl2,'TableGPSAll_8.txt','Delimiter',',')
writetable(Tabl3,'TableGPSAll_9.txt','Delimiter',',')
writetable(Tabl4,'TableGPSAll_10.txt','Delimiter',',')



% ---------Add Elevation data from Movebank: ASTER Global DEM A-----------
% for the updated file 
load 1HzFirstAnalysis_NewTimeString
TableGPSAll=[TableGPSAll(:,1:16),TableGPSAll(:,19)];
% read the data
load 1HzFirstAnalysis_Spring
Ale1_1=readtable('WithEnvData\Anotation Elevation\TableGPSAll_1_EnvEle.csv');
Ale1_2=readtable('WithEnvData\Anotation Elevation\TableGPSAll_2_EnvEle.csv');
Ale1_3=readtable('WithEnvData\Anotation Elevation\TableGPSAll_3_EnvEle.csv');
Ale1_4=readtable('WithEnvData\Anotation Elevation\TableGPSAll_4_EnvEle.csv');
Ale1_5=readtable('WithEnvData\Anotation Elevation\TableGPSAll_5_EnvEle.csv');
Ale1_6=readtable('WithEnvData\Anotation Elevation\TableGPSAll_6_EnvEle.csv');

AleAboveTerrain=[Ale1_1.ASTERASTGTM2Elevation;Ale1_2.ASTERASTGTM2Elevation;...
    Ale1_3.ASTERASTGTM2Elevation;Ale1_4.ASTERASTGTM2Elevation;Ale1_5.ASTERASTGTM2Elevation;...
   Ale1_6.ASTERASTGTM2Elevation];

ELE=TableGPSAllNew.InterpolatedElevation-AleAboveTerrain;
ELE(ELE<0)=0;
TableGPSAllNew.AleAboveTerrain=ELE;
TableGPSAllNew.TerrainHeight=AleAboveTerrain;
save('1HzGPS_and_ACC_final','TableGPSAllNew','SensorsAll','StatsTCleaned')
save('StatsTCleaned','StatsTCleaned')



% ---------Add Elevation data from Movebank: ASTER Global DEM B-----------
% for the updated file 
load 1HzGPS_and_ACC_2update
TableGPSAllNew_2add.running_Flap_rate_4sec(:,1)=nan;
% read the data
Ale1_1=readtable('WithEnvData\Anotation Elevation\TableGPSAll_7_EnvEle.csv');
Ale1_2=readtable('WithEnvData\Anotation Elevation\TableGPSAll_8_EnvEle.csv');
Ale1_3=readtable('WithEnvData\Anotation Elevation\TableGPSAll_9_EnvEle.csv');
Ale1_4=readtable('WithEnvData\Anotation Elevation\TableGPSAll_10_EnvEle.csv');

AleAboveTerrain=[Ale1_1.ASTERASTGTM2Elevation;Ale1_2.ASTERASTGTM2Elevation;...
    Ale1_3.ASTERASTGTM2Elevation;Ale1_4.ASTERASTGTM2Elevation];

ELE=TableGPSAllNew_2add.InterpolatedElevation-AleAboveTerrain;
ELE(ELE<0)=0;
TableGPSAllNew_2add.AleAboveTerrain=ELE;
TableGPSAllNew_2add.TerrainHeight=AleAboveTerrain;
save('1HzGPS_and_ACC_2update','TableGPSAllNew_2add','SensorsAll_2add','StatsTCleaned_2add')
save('StatsTCleaned_2update','StatsTCleaned_2add')