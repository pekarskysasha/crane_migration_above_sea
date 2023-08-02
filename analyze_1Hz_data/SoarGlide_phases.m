function [Answer,soar_glideWork]=SoarGlide_phases(soar_glide,TableGPSCont,indev)
%% --(1) Relating a glide to a cisling phase
%-----(a) only the first found glide
%-----(b) only if it starts not more them 60 sec after soar end
[~,RelevantIndexes]=unique(soar_glide(:,1)); % only the first glide
soar_glideFirst=soar_glide(RelevantIndexes,:);
IndexRelInterval=(soar_glideFirst(:,3)-soar_glideFirst(:,2))<61; % relevamt interval
soar_glideWork=soar_glideFirst(IndexRelInterval,:);
if isempty(soar_glideWork)
    Answer=[];
    soar_glideWork=[];
    return;
end
%% --(2) glide and soar stats
s=size(soar_glideWork,1);
soart=soar_glideWork(:,2)-soar_glideWork(:,1); % soar duration
glidet=soar_glideWork(:,4)-soar_glideWork(:,3); % glide duration
climbrate=(TableGPSCont.InterpolatedElevation(soar_glideWork(:,2))-TableGPSCont.InterpolatedElevation(soar_glideWork(:,1)))./soart(end-s+1:end);
ExitAlt=TableGPSCont.AleAboveTerrain(soar_glideWork(:,2)); % soar exit altitude

%% --(3) Flight parameters
%--------(a) sink speed during glide m/sec
vz=(TableGPSCont.InterpolatedElevation(soar_glideWork(:,4))-TableGPSCont.InterpolatedElevation(soar_glideWork(:,3)))./glidet(end-s+1:end);
%--------(b) ground speed during glide m/sec
vg=hypot(TableGPSCont.Nir_X(soar_glideWork(:,4))-TableGPSCont.Nir_X(soar_glideWork(:,3)),TableGPSCont.Nir_Y(soar_glideWork(:,4))-TableGPSCont.Nir_Y(soar_glideWork(:,3)))./glidet(end-s+1:end);
%--------(c) distance traveled when gliding
gl_dist_m=round(hypot(TableGPSCont.Nir_X(soar_glideWork(:,4))-TableGPSCont.Nir_X(soar_glideWork(:,3)),TableGPSCont.Nir_Y(soar_glideWork(:,4))-TableGPSCont.Nir_Y(soar_glideWork(:,3))));
%--------(d) air speed and tail wind during glide
ang=atan2(TableGPSCont.Nir_Y(soar_glideWork(:,4))-TableGPSCont.Nir_Y(soar_glideWork(:,3)),TableGPSCont.Nir_X(soar_glideWork(:,4))-TableGPSCont.Nir_X(soar_glideWork(:,3)));
%-- (d.1)--- Wind from ECMWF-----------------
wu=[];
wv=[];
for j=1:s
    wu(j,1)=mean(TableGPSCont.u925(soar_glideWork(j,3):soar_glideWork(j,4)));
    wv(j,1)=mean(TableGPSCont.v925(soar_glideWork(j,3):soar_glideWork(j,4)));
end
% atan2(Y,X)
% V => Velocity of the north-south (meridoinal) component of wind. Positive values indicate south to north flow (m/s) => Y
% U => Velocity of the east-west (zonal) component of wind. Positive values indicate west to east flow (m/s) => X
wang=atan2(wv,wu);
wvel=hypot(wu,wv);
% tail wind
tw=cos(wang-ang).*wvel;
% side wind
sw=abs(sin(wang-ang)).*wvel;
% air speed
va=hypot(vg(end-s+1:end)-tw(end-s+1:end),sw(end-s+1:end));

%-- (c.2)--- Wind from thermal drift-----------
Indexes=soar_glideWork(:,1:2);
[wang_calc,wvel_calc]=WindCalcultedThermals(Indexes,TableGPSCont);
%--- if the thermal us shorter then 72 sec, wind not calculated, fill from average
%--check untervals to the existing thermal winds, and if possible (15 min intervals) fill in
IndNan=find(isnan(wang_calc));
RunIndex=[1:1:length(wang_calc)]';
wang_calc_filled=wang_calc;
wvel_calc_filled=wvel_calc;
for f=1:length(IndNan)
    ind4Average=[];
    if IndNan(f)==1 % if first thermal, take next existing
        indNext=find(~isnan(wang_calc) & RunIndex>=IndNan(f)+1,1,'first');
        ind4Average=[indNext,Indexes(indNext,1)-Indexes(IndNan(f),2)];
    elseif IndNan(f)==length(wang_calc) % if first thermal, thake previous existing
        indPrev=find(~isnan(wang_calc) & RunIndex<=IndNan(f)-1,1,'last');
        ind4Average=[indPrev,Indexes(IndNan(f),1)-Indexes(indPrev,2)];
    else %thake previous and next
        indNext=find(~isnan(wang_calc) & RunIndex>=IndNan(f)+1,1,'first');
        ind4Average=[indNext,Indexes(indNext,1)-Indexes(IndNan(f),2)];
        indPrev=find(~isnan(wang_calc) & RunIndex<=IndNan(f)-1,1,'last');
        ind4Average=[ind4Average; [indPrev,Indexes(IndNan(f),1)-Indexes(indPrev,2)]];
    end
    % fill the mean wind values if the theramls are no more then 15 min apart
    if ~isempty(ind4Average)
        wang_calc_filled(IndNan(f))=nanmean(wang_calc(ind4Average(ind4Average(:,2)<901),1));
        wvel_calc_filled(IndNan(f))=nanmean(wvel_calc(ind4Average(ind4Average(:,2)<901),1));
    end
end
wang_calc=wang_calc_filled;
wvel_calc=wvel_calc_filled;
% tail wind
tw_calc=cos(wang_calc-ang).*wvel_calc;
% side wind
sw_calc=abs(sin(wang_calc-ang)).*wvel_calc;
% air speed
va_calc=hypot(vg(end-s+1:end)-tw_calc(end-s+1:end),sw_calc(end-s+1:end));


%% --(3) Atmospheric conditions and flap rate
sfr=nan(s,1);
gfr=nan(s,1);
sfp=nan(s,1);
gfp=nan(s,1);
blh=nan(s,1);
sl=nan(s,1);
day_time=nan(s,1);
for j=1:s
    %--------(a) mean flap rate during climb and during glide
    sfr(j)=nanmean(TableGPSCont.running_Flap_rate_4sec(soar_glideWork(j,1):soar_glideWork(j,2)));
    gfr(j)=nanmean(TableGPSCont.running_Flap_rate_4sec(soar_glideWork(j,3):soar_glideWork(j,4)));
    %--------(b) proportion of seconds with flapping
    PlapNum=round(TableGPSCont.running_Flap_rate_4sec(soar_glideWork(j,1):soar_glideWork(j,2)));
    sfp(j)=sum(PlapNum>0)/length(PlapNum);
    PlapNum=round(TableGPSCont.running_Flap_rate_4sec(soar_glideWork(j,3):soar_glideWork(j,4)));
    gfp(j)=sum(PlapNum>0)/length(PlapNum);
    %--------(c) boundry layer hight
    blh(j)=mean(TableGPSCont.blh(soar_glideWork(j,3):soar_glideWork(j,4)));
    %--------(d) sea or land (if both leave nan)
    UniqueSeaLand=unique(TableGPSCont.OverSeaOrLand(soar_glideWork(j,1):soar_glideWork(j,4)));
    UniqueSeaLand=UniqueSeaLand(~isnan(UniqueSeaLand));
    if length(UniqueSeaLand)==1
        sl(j)=UniqueSeaLand;
    end
    %--------(e) Time of day (if any point during night, define as night)
    if sum(TableGPSCont.day_light(soar_glideWork(j,1):soar_glideWork(j,4))==0)>0
        day_time(j)=0;
    else
        day_time(j)=1;
    end  
end
%% --(4) create the final table
DateTimeSrart=TableGPSCont.TimeCont(soar_glideWork(:,1));
Individual=ones(length(soart),1)*indev;
lon=TableGPSCont.Interpolated_Lon(soar_glideWork(:,1));
lat=TableGPSCont.Interpolated_Lat(soar_glideWork(:,1));
if size(soar_glideWork,2)==5
    ThermalType=soar_glideWork(:,5);
else
    ThermalType=ones(size(soar_glideWork,1),1);
end
Answer=table(Individual,DateTimeSrart,lon,lat,soart,climbrate,sfr,sfp,ThermalType,ExitAlt,glidet,gl_dist_m,vz,va,va_calc,...
    tw,tw_calc,sw,sw_calc,gfr,gfp,blh,sl,day_time,...
    'VariableNames',{'individual','datetime_start','lon_start','lat_start','time_soaring_sec','climb_rate_m_sec',...
    'falp_rate_climb','falp_prop_climb','type_thermal','exit_alt_above_terr','time_gliding_sec','distance_gliding_m',...
    'sink_speed_m_sec','air_speed_m_sec','air_speed_calc_m_sec','tail_wind','tail_wind_calc','side_wind','side_wind_calc',...
    'falp_rate_glide','falp_prop_glide','mean_blh','sea_land','day_time'});
end