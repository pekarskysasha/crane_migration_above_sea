% clear all
addpath('ThermalAnalysis')
addpath('ThermalAnalysis/Working folder/FinalTracks')
load ToEraze
load StatsTCleaned % meta data of the clened parts
Indiv=unique(StatsTCleaned.Individual);
% for looking for full trips
load MigrationAnalysis_Stage3_TripLevel
% for age and sex data
load WorkingMetaData
% for knowing the ditance to coast
[Hsea, Asea] = shaperead('../../LandCover\Seas\World_Seas_IHO_v3','UseGeoCoords', true);
MS_EB=table(Hsea(44).Lon', Hsea(44).Lat','VariableNames',{'Lon','Lat'}); % Mediterranean Sea - Eastern Basin
BS=table(Hsea(90).Lon', Hsea(90).Lat','VariableNames',{'Lon','Lat'}); %Black sea
%% the directory with the data per crane
directory='ThermalAnalysis/Working folder/FinalTracks';
FilesIn=struct2table(dir(directory));
IndevidualsExistTemp=FilesIn.name(FilesIn.bytes>0);
IndevidualsExist=[];
Dist2coastMed=[];
Dist2coastBlack=[];
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
CiclingThemalOnly=[];
ID=0;
for i=1:length(Indiv)
    Animal_ID=Indiv(i);
    indArchaive=find(IndevidualsExist==Animal_ID);
    if isempty(indArchaive)
        continue;
    end
    load([directory,'\1HzGPS_and_ACC_annotated_',num2str(IndevidualsExist(indArchaive)),'.mat']);
    %TableGPSAllNew=eval(['TableGPSAllNew_',mat2str(Animal_ID)]);
    AllInfo=StatsTCleaned(StatsTCleaned.Individual==Indiv(i),:);
    for h=1:height(AllInfo)  
        % (I)-- find age and sex --------------------------------
            %---find in metadata
            ind=MetaData.Animal_ID==Animal_ID;
            % -- Crane status----------------------------------------------------
            if MetaData.Is_Breeding_Adult(ind)==1 % breeding adult
                Status=1;
                StatusS={'breeding adult'};
            else
                % date of hatch
                if strcmp(MetaData.Trapping_Age(ind),'Adult')
                    if month(MetaData.Date_deployed_on(ind))<3 % tagged spring
                        YearHatch=year(MetaData.Date_deployed_on(ind))-3;
                    else % tagged fall
                        YearHatch=year(MetaData.Date_deployed_on(ind))-2;
                    end                   
                elseif strcmp(MetaData.Trapping_Age(ind),'Juvenile')
                    if month(MetaData.Date_deployed_on(ind))<3 % tagged spring
                        YearHatch=year(MetaData.Date_deployed_on(ind))-1;
                    else % tagged fall
                        YearHatch=year(MetaData.Date_deployed_on(ind));
                    end
                elseif strcmp(MetaData.Trapping_Age(ind),'Subadult')
                    if month(MetaData.Date_deployed_on(ind))<3 % tagged spring
                        YearHatch=year(MetaData.Date_deployed_on(ind))-2;
                    else % tagged fall
                        YearHatch=year(MetaData.Date_deployed_on(ind))-1;
                    end
                end
                
                AgeNum=AllInfo.Year(h)-YearHatch;
                if AgeNum>1 & strcmp(AllInfo.Season(h),'Fall')
                    Status=2;
                    StatusS={'Unknown adult'};
                elseif AgeNum>2 & strcmp(AllInfo.Season(h),'Spring')
                    Status=2;
                    StatusS={'Unknown adult'};
                elseif AgeNum==0 & strcmp(AllInfo.Season(h),'Fall')
                    Status=4;
                    StatusS={'Juvenile'}; 
                elseif AgeNum==1 & strcmp(AllInfo.Season(h),'Spring')
                    Status=4;
                    StatusS={'Juvenile'};   
                elseif AgeNum==1 & strcmp(AllInfo.Season(h),'Fall')
                    Status=3;
                    StatusS={'Subadult'};  
                 elseif AgeNum==2 & strcmp(AllInfo.Season(h),'Spring')
                    Status=3;
                    StatusS={'Subadult'};    
                end
            end
            age=Status;
            age_w=StatusS;
            % -- Crane sex ----------------------------------------------------
            sex=MetaData.Sex(ind);
        
        % (II)-- calculate thermal cicling and gliding --------------------------------
        SG=[];
        ThermalTableNirF=[];
        disp(['working on individual: ',num2str(Animal_ID),' section starts on ',char(AllInfo.TimeStartV(h))])
        TableGPSCont=TableGPSAllNew(TableGPSAllNew.UniqueSectionCounter==AllInfo.UniqueSectionCounter(h),:);
        %________________________________________________________________________________________________
        
        %% -(1)---Find thermals and Gliding using Nir's code
        %----(a) find indexes of circling/gliding
        %----(b) create information table for cicling soaring phase
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        InterpolatedElevation=TableGPSCont.InterpolatedElevation;
        [soar_phase,soar_glide]=FindSoarAndGlide(TableGPSCont.Nir_X,TableGPSCont.Nir_Y,...
            InterpolatedElevation,datenum(TableGPSCont.TimeCont),TableGPSCont.mag_x);
        %         circ=FindSoaringCircling1(TableGPSCont.Nir_X,TableGPSCont.Nir_Y,InterpolatedElevation);
        %         circ=FindSoaringCirclingOrNot1(TableGPSCont.Nir_X,TableGPSCont.Nir_Y,InterpolatedElevation);
        if ~isempty(soar_phase)
            %             [soar_phase,soar_glide]=FindGlide1_old(TableGPSCont.Nir_X,TableGPSCont.Nir_Y,InterpolatedElevation,datenum(TableGPSCont.TimeCont),circ);
            %             [soar_phase,soar_glide]=FindGlide1(TableGPSCont.Nir_X,TableGPSCont.Nir_Y,InterpolatedElevation,datenum(TableGPSCont.TimeCont),circ);
            [ThermalTableNir]=OrganiseThermalTable(soar_phase,TableGPSCont.TimeCont,InterpolatedElevation,Animal_ID);
        else
            soar_glide=[];
            ThermalTableNir=[];
        end
        if ~isempty(ThermalTableNir)
            ThermalTableNirF=AddInfoToThermalTable(ThermalTableNir,TableGPSCont);
            [wang_calc,wvel_calc]=WindCalcultedThermals(soar_phase,TableGPSCont);
            ThermalTableNirF=[ThermalTableNirF,table(wang_calc,wvel_calc)];
      %% ==== remove lines that were manualy found to be a mistake =====
            [k,ID_Remove]=intersect(ThermalTableNirF(:,[1:3,6:9]),ToEraze(:,[1:3,6:9]));
            if ~isempty(ID_Remove)
                ThermalTableNirF(ID_Remove,:)=[];
                ThermalTableNir(ID_Remove,:)=[];
                for zz=1:length(ID_Remove)
                    IDRemove_soar_phase=find(soar_phase(:,1)==k.IndexStart(zz) & soar_phase(:,2)==k.IndexEnd(zz));
                    soar_phase(IDRemove_soar_phase,:)=[];
                    if ~isempty(soar_glide)
                        IDRemove_soar_glide=find(soar_glide(:,1)==k.IndexStart(zz) & soar_glide(:,2)==k.IndexEnd(zz));
                        soar_phase(IDRemove_soar_glide,:)=[];
                    end
                end
            end
      %%================================================================    
        end
        CiclingThemalOnly=[CiclingThemalOnly; ThermalTableNirF];
        %________________________________________________________________________________________________
%        if ~isempty(ThermalTableNirF) & sum(ThermalTableNirF.AverageTimePerTurn>150)>0
%            ThermalTableNirF;
%        end
%       
        %% -(2)---Alocate gliding phase to thermal cicling--------------------------------
        %-----(a) only the first found glide
        %-----(b) only if it starts not more them 60 sec after soar ends
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isempty(soar_glide)
            [SG,UsedIndexes4Plot]=Soar_glide_phases(soar_glide,TableGPSCont,Animal_ID);
            SG.age(:,1)=age;
            SG.age_w(:,1)=age_w;
            SG.sex(:,1)=sex;
            SG_all=[SG_all; SG];
            
            %             if sum(SG.exit_alt_above_terr>3000)>0 % | sum(SG.distance_gliding_m./SG.time_soaring_sec>50)>0
            %                 SG;
            %             end
            %
            %             X = rand;
            %             if X<0.15 & sum(SG.sea_land~=0)>0
            %                 X;
            %             end
            %============= look for sections to plot ===============================
            %-- (1) more then 5 turns for a soar with a glide
            %-- (2) at least 3 soar glides
            %-- (3) at least 6 soars
            %-- (4) total duration longer than 40
            %-- (5) look for ditance to coast for sea soring
            if ~isempty(intersect(soar_phase(soar_phase(:,4)>5,1),soar_glide(:,1))) &...
                    size(soar_glide,1)>2 & size(soar_phase,1)>5 & AllInfo.TimeOfBurstT(h)>=40
                %-- (5) the trip is full
                indInfoTrip=find(MigrationInfo.animal_ID==Animal_ID ...
                    & datenum(MigrationInfo.start_mig)<=floor(datenum(TableGPSCont.TimeCont(1))) &...
                    datenum(MigrationInfo.end_mig)>=floor(datenum(TableGPSCont.TimeCont(1))));
                
                if isempty(indInfoTrip) | MigrationInfo.complete_track(indInfoTrip)==0 |  MigrationInfo.N_missing_days(indInfoTrip)> 4
                    continue;
                end
                
                if TableGPSCont.Latitude(1)>49 %(a) north
                    AllInfo.UniqueSectionCounter(h);
                elseif sum(SG.sea_land==1)>0 %(b) black sea
                    AllInfo.UniqueSectionCounter(h);                    
                    % distance to coast
                    for s=1:height(SG)
                        Dist2coast=deg2km(p_poly_dist(SG.lat_start(s),SG.lon_start(s),BS.Lat,BS.Lon));
                        if Dist2coast<0
                            if Dist2coast < -120
                                Dist2coast
                            end
                            Dist2coastBlack=[Dist2coastBlack; Dist2coast];
                        end
                    end
                elseif sum(SG.sea_land==2)>0 %(c) med sea
                    AllInfo.UniqueSectionCounter(h);
                    % distance to coast
                    for s=1:height(SG)
                        Dist2coast=deg2km(p_poly_dist(SG.lat_start(s),SG.lon_start(s),MS_EB.Lat,MS_EB.Lon));
                        if Dist2coast<0
                            if Dist2coast < -55
                                Dist2coast
                            end
                            Dist2coastMed=[Dist2coastMed; Dist2coast];
                        end
                    end
                elseif TableGPSCont.Latitude(1)<25 %(d) south
                    AllInfo.UniqueSectionCounter(h);
                end
                %PlotSoaringGliding_4NIR(AllInfo(h,:),TableGPSCont,SG,ThermalTableNirF,AllInfo.UniqueSectionCounter(h))   
            end
            
        end
        if sum(soar_phase(:,3)==2)>0
            if sum(soar_phase(:,3)==2 & soar_phase(:,4)==2)>0
                soar_phase;
            elseif sum(soar_phase(:,3)==2 & soar_phase(:,4)>2)>0
                soar_phase;
            end
        end

        %% -(3)--- Make stats for thermals only (No glide)
        %---Using the ThermalTable for that
        % --------Analyse thermals only for climb rate and other parameters ---------
        %-----(a) sections length > 10 minutes
        %-----(b) the data devided into 10 minute sections
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [TheramlSatats]=equal_interval_thermal(Animal_ID,ThermalTableNir,TableGPSCont);
        if ~isempty(TheramlSatats)
            TheramlSatats.age(:,1)=age;
            TheramlSatats.age_w(:,1)=age_w;
            TheramlSatats.sex(:,1)=sex;
        end
        TheramlStatsAll=[TheramlStatsAll; TheramlSatats];
        if ~isempty(TheramlSatats) & sum(TheramlSatats.time_sice_sunset_h>1 & TheramlSatats.percent_time_in_thermals>0)>0
           TheramlSatats; 
        end
        %% -(4)---  Analyse flap rate at different stages------------------------------
        %-the analysis on all the thermals (from thermal table and not from soar-glide) and also random sections of 5 minutes from 1Hz parts with no soaring at all
        %--CircleORNotColum:
        %------[1] is thermal
        %------[0] is not thermal
        %------[2] No thermals, flapping flight
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isempty(ThermalTableNir)
            [FlapT]=AnalyseFlapRate(ThermalTableNir,TableGPSCont.Interpolated_Lat,TableGPSCont.Interpolated_Lon,...
                TableGPSCont.TimeCont,TableGPSCont.AleAboveTerrain,TableGPSCont.running_Flap_rate_4sec,TableGPSCont.OverSeaOrLand);
            if ~isempty(FlapT)
                % associate each circling event with its corresponding between-cicling (proxy for glide) section
                FlapT.runing_ID=FlapT.runing_ID+ID;
                ID=FlapT.runing_ID(end);
                FlapAll=[FlapAll; FlapT];
            end
        else
            %-- Create random indexes of 5 minute length for making flap rate where there is no soaring
            IndexRandom=[];
            indexEnd=0;
            while indexEnd<height(TableGPSCont)-321
                IndexRandom=[IndexRandom; [indexEnd+randi(20, [1,1]),indexEnd+randi(20, [1,1])+300]];
                indexEnd=indexEnd+randi(20, [1,1])+300;
            end
            IndexRandom=table(IndexRandom(:,1),IndexRandom(:,2),'VariableNames',{'IndexStart','IndexEnd'});
            [FlapNotSoar]=AnalyseFlapRateNoSoar(IndexRandom,TableGPSCont.Interpolated_Lat,TableGPSCont.Interpolated_Lon,...
                TableGPSCont.TimeCont,TableGPSCont.AleAboveTerrain,TableGPSCont.running_Flap_rate_4sec,TableGPSCont.OverSeaOrLand);
            if ~isempty(FlapNotSoar)
                FlapNotSoar.runing_ID=FlapNotSoar.runing_ID+ID;
                ID=FlapNotSoar.runing_ID(end);
                FlapAll=[FlapAll; FlapNotSoar];
            end
        end
        
         if ~isempty(TheramlSatats) 
             IndStrangeThermal=find(TheramlSatats.OverSea>0 & TheramlSatats.percent_time_in_thermals>0);
             if ~isempty(IndStrangeThermal) %&& sum(TheramlSatats.Mean_DeltaT(IndStrangeThermal)<0)>0
                 ThermalTableNirF;
             end
         end
        %% -(5)---  Various plotting function for data exploration -------------------------
        %plotTracks(TableGPSCont,soar_glide,UsedIndexes4Plot,soar_phase)
        %plotCirc(circ,TableGPSCont.Nir_X,TableGPSCont.Nir_Y)
    end
end
%% ---- Analyse the soaring with gliding --------------------------------
%------Plot glide polar -------------------------------------------
figure
Vair=9:1:30; % True air speed for which sink speed is given
Vsink=-[0.668,0.6627,0.6774,0.7038,0.7409,0.7881,0.845,0.9112,0.9865,1.071,1.164,1.265,1.375,1.493,1.62,1.755,1.898,2.05,2.211,2.381,2.563,2.756];
scatter(SG_all.air_speed_m_sec,SG_all.sink_speed_m_sec,20,'MarkerEdgeColor','k','MarkerFaceColor',[162, 166, 173]/255)
hold on
plot(Vair,Vsink,'-r')
ylabel('Sink rate (m sec^-^1)')
xlabel('Air speed (m sec^-^1)')
ax1 = gca;
ax1.FontSize=16;

%-------------- For V-opt, given real climb rate -------------------------------------------
OrgClr=0.5:0.5:3; % Climb in thermals,descrete values foe which Vopt is given
OrgVopt=[17.9,21,23.7,26,28,29.5]; % Optimum inter-thermal speed (given as descrete values)
Vopt=interp1(OrgClr,OrgVopt,SG_all.climb_rate_m_sec,'spline'); % Interpolated Vopt in better resolution
Vbg=14.4;

figure
scatter(SG_all.climb_rate_m_sec,SG_all.sink_speed_m_sec,20,'MarkerEdgeColor','k','MarkerFaceColor',[162, 166, 173]/255)
%-------------- RAFI ------------------------------------------------------------------------
RAFI=(Vopt-SG_all.air_speed_m_sec)./(Vopt-ones(length(Vopt),1)*Vbg);
RAFI_calc=(Vopt-SG_all.air_speed_calc_m_sec)./(Vopt-ones(length(Vopt),1)*Vbg);
%-------------- Save -----------------------------------------------
TEA=SG_all.exit_alt_above_terr;
SoarGlide=SG_all;
SoarGlide.Vopt=Vopt;
SoarGlide.RAFI=RAFI;
SoarGlide.RAFIcalc=RAFI_calc;
SoarGlide.Date=datetime(SoarGlide.datetime_start,'format','dd-MMM-y');
writetable(CiclingThemalOnly,'SoarGlide_new_FlapRationNew.csv')
writetable(TheramlStatsAll,'TheramlStats_new.csv')
writetable(FlapAll,'FlapAll_new.csv')
writetable(CiclingThemalOnly,'ThermalCiclingOnly_new.csv')
save('Soar_new','SoarGlide','TheramlStatsAll','FlapAll','CiclingThemalOnly')

%% =========== Make a ThermalTable from Nir's soar_phase indexes ==========
function ThermalTable=OrganiseThermalTable(soar_phase,T,InterpolatedElevation,Animal_ID)
%-Take Nir's soar_phase inedexes and calculates the stats
Times=[];
if ~isempty(soar_phase)
    for i=1:length(soar_phase(:,1))
        IndexStart=soar_phase(i,1);
        IndexEnd=soar_phase(i,2);
        if size(soar_phase,2)==3
            NumberofTurns=soar_phase(i,3);
            ThermalType=1;
        else
            NumberofTurns=soar_phase(i,4);
            ThermalType=soar_phase(i,3);
        end
        TimeInThermal=datenum((T(IndexEnd)-T(IndexStart)))*86400; % sec
        ClimbRate=(InterpolatedElevation(IndexEnd)-InterpolatedElevation(IndexStart))/TimeInThermal;
        AverageTimePerTurn=(InterpolatedElevation(IndexEnd)-InterpolatedElevation(IndexStart))/NumberofTurns;
        MeanElevationGainPerTurn=(InterpolatedElevation(IndexEnd)-InterpolatedElevation(IndexStart))/NumberofTurns;
        
        Idexe(i,:)=[Animal_ID,IndexStart,IndexEnd,TimeInThermal,NumberofTurns,AverageTimePerTurn,MeanElevationGainPerTurn,ClimbRate,ThermalType];
        Times=[Times; [T(IndexStart),T(IndexEnd)]];
    end
    ThermalTable=table(Idexe(:,1),Idexe(:,2),Idexe(:,3),Times(:,1),Times(:,2),Idexe(:,4),Idexe(:,5),Idexe(:,6),Idexe(:,7),Idexe(:,8),Idexe(:,9),...
        'VariableNames',{'Animal_ID','IndexStart','IndexEnd','TimeStart','TimeEnd','TimeInThermal',...
        'NumberOfTurns','AverageTimePerTurn','MeanElevationGainPerTurn','ClimbRate', 'ThermalType'});
else
    ThermalTable=table([]);
end
end

%% =================== add info to the thermal table for analysis ===============================
function [ThermalTableFinal]=AddInfoToThermalTable(ThermalTable,TableGPS)

for i=1:height(ThermalTable)
    max_elvation_above_ground(i,1)=nanmax(TableGPS.AleAboveTerrain(ThermalTable.IndexStart(i):ThermalTable.IndexEnd(i)));
    mean_blh(i,1)=nanmean(TableGPS.blh(ThermalTable.IndexStart(i):ThermalTable.IndexEnd(i)));
    mean_Flap_rate(i,1)=nanmean(TableGPS.running_Flap_rate_4sec(ThermalTable.IndexStart(i):ThermalTable.IndexEnd(i)));
    UniqueSeaLand=unique(TableGPS.OverSeaOrLand(ThermalTable.IndexStart(i):ThermalTable.IndexEnd(i)));
    UniqueSeaLand=UniqueSeaLand(~isnan(UniqueSeaLand));
    if length(UniqueSeaLand)==1
        OverSeaOrLand(i,1)=UniqueSeaLand;
    else
        OverSeaOrLand(i,1)=nan;
    end
    OverSeaOrLand(i,1)=nanmedian(TableGPS.OverSeaOrLand(ThermalTable.IndexStart(i):ThermalTable.IndexEnd(i)));
    if ~isnan(TableGPS.Interpolated_Lat(ThermalTable.IndexStart(i)))
        lat_start(i,1)=TableGPS.Interpolated_Lat(ThermalTable.IndexStart(i));
        lon_start(i,1)=TableGPS.Interpolated_Lon(ThermalTable.IndexStart(i));
    else
        lat_start(i,1)=TableGPS.Interpolated_Lat(ThermalTable.IndexEnd(i));
        lon_start(i,1)=TableGPS.Interpolated_Lon(ThermalTable.IndexEnd(i));
    end
    % if any point during night, define as night
    if sum(TableGPS.day_light(ThermalTable.IndexStart(i):ThermalTable.IndexEnd(i))==0)>0
        day_time(i,1)=0;
    else
        day_time(i,1)=1;
    end
    %-- (a) flight parameters
    wu=mean(TableGPS.u925(ThermalTable.IndexStart(i):ThermalTable.IndexEnd(i)));
    wv=mean(TableGPS.v925(ThermalTable.IndexStart(i):ThermalTable.IndexEnd(i)));
    
    % atan2(Y,X)
    % V => Velocity of the north-south (meridoinal) component of wind. Positive values indicate south to north flow (m/s) => Y
    % U => Velocity of the east-west (zonal) component of wind. Positive values indicate west to east flow (m/s) => X
    wang(i,1)=atan2(wv,wu);
    wvel(i,1)=hypot(wu,wv);
end
ThermalTableFinal=[ThermalTable, table(lat_start,lon_start,max_elvation_above_ground,mean_blh,wang,wvel,OverSeaOrLand,mean_Flap_rate,day_time)];
end

%% ===========Plot function: Elevation + Magnetometer =====================
function plotTracks(TableGPS,soar_glide,SaorGlideTable,soar_phase)
Ele=TableGPS.InterpolatedElevation;
InterpolatedElevation=Ele;
X=[1:1:length(TableGPS.TimeCont)];
TH=TableGPS.TerrainHeight;
if sum(~isnan(TH))==0
    TH=zeros(length(TH),1);
end
IntrTH = interp1(X(~isnan(TH)),TH(~isnan(TH)),X);
figure;

%% -------------- Raw analysis----------------------------------------
subplot(2,1,1)
hold on
ClimbGlideIndex=zeros(height(TableGPS),1);
ForPlot=Ele;
for e=1:length(soar_phase(:,1))
    if size(soar_phase,2)==4
        if soar_phase(e,3)==1
            ClimbGlideIndex(soar_phase(e,1):soar_phase(e,2))=1;
        elseif soar_phase(e,3)==2
            ClimbGlideIndex(soar_phase(e,1):soar_phase(e,2))=2;
        else
            ClimbGlideIndex(soar_phase(e,1):soar_phase(e,2))=3;
        end
    else
        ClimbGlideIndex(soar_phase(e,1):soar_phase(e,2))=1;
    end
end

for p=1:length(soar_glide(:,1))
    S1=soar_glide(p,1);
    E1=soar_glide(p,2);
    S2=soar_glide(p,3);
    E2=soar_glide(p,4);
    
    yyaxis left
    plot(X(S1:E1),InterpolatedElevation(S1:E1),'-','Color',[249, 117, 216]/255,'LineWidth',3);
    ForPlot(S1:E1)=nan;
    if ~isnan(S2)
        plot(X(S2:E2),InterpolatedElevation(S2:E2),'-','Color',[73, 145, 193]/255,'LineWidth',3);
        ForPlot(S2:E2)=nan;
    end
end
plot(X,ForPlot,'-','Color','k','LineWidth',2);
area(X,IntrTH,'FaceColor',[188, 146, 101]/255)
ax1 = gca;
ax1.YLim=[0 max(InterpolatedElevation)+100];
ylabel('Elevation above sea level(m)')

yyaxis right
Z1=TableGPS.mag_x;
Z1(ClimbGlideIndex==0 | ClimbGlideIndex==2 | ClimbGlideIndex==3)=nan; % 1
Z2=TableGPS.mag_x;
Z2(ClimbGlideIndex==1 | ClimbGlideIndex==2 | ClimbGlideIndex==3)=nan; % 0
Z3=TableGPS.mag_x;
Z3(ClimbGlideIndex==1 | ClimbGlideIndex==3 | ClimbGlideIndex==0)=nan; % 2
Z4=TableGPS.mag_x;
Z4(ClimbGlideIndex==1 | ClimbGlideIndex==2 | ClimbGlideIndex==0)=nan; % 3

P1=plot(X,Z1,'-','Color',[232, 99, 46]/255,'LineWidth',0.3);
hold on
P1=plot(X,Z3,'-','Color',[88, 128, 214]/255,'LineWidth',0.3);
P1=plot(X,Z4,'-','Color',[159, 73, 209]/255,'LineWidth',0.3);
P2=plot(X,Z2,'-','Color',[90, 87, 92]/255,'LineWidth',0.3);
P1.Color(4) = 0.5;
P2.Color(4) = 0.5;
%% --------- Genral Title--------------------------------------------
title('Raw soar and glide on elevation, mag. color by soring only')

datacursormode on
dcm = datacursormode(gcf);
set(dcm,'UpdateFcn',@myupdatefcn)
%% -------------- Soar glide after final analysis----------------------------------------
subplot(2,1,2)
hold on
ClimbGlideIndex=zeros(height(TableGPS),1);
ForPlot=Ele;
for e=1:length(SaorGlideTable(:,1))
    ClimbGlideIndex(SaorGlideTable(e,1):SaorGlideTable(e,2))=1;
end

for p=1:length(SaorGlideTable(:,1))
    S1=SaorGlideTable(p,1);
    E1=SaorGlideTable(p,2);
    S2=SaorGlideTable(p,3);
    E2=SaorGlideTable(p,4);
    
    yyaxis left
    plot(X(S1:E1),InterpolatedElevation(S1:E1),'-','Color',[249, 117, 216]/255,'LineWidth',3);
    ForPlot(S1:E1)=nan;
    if ~isnan(S2)
        plot(X(S2:E2),InterpolatedElevation(S2:E2),'-','Color',[73, 145, 193]/255,'LineWidth',3);
        ForPlot(S2:E2)=nan;
    end
end
plot(X,ForPlot,'-','Color','k','LineWidth',2);
area(X,IntrTH,'FaceColor',[188, 146, 101]/255)
ax1 = gca;
ax1.YLim=[0 max(InterpolatedElevation)+100];
ylabel('Elevation above sea level(m)')

yyaxis right
Z1=TableGPS.mag_x;
Z1(ClimbGlideIndex==0)=nan;
Z2=TableGPS.mag_x;
Z2(ClimbGlideIndex==1)=nan;
P1=plot(X,Z1,'-','Color','r','LineWidth',0.3);
hold on
P2=plot(X,Z2,'-','Color','b','LineWidth',0.3);
P1.Color(4) = 0.5;
P2.Color(4) = 0.5;
%% ---------Title--------------------------------------------
title('final (analysed) soaring-gliding on elevation, mag. color soar phase of final soaring-gliding')

datacursormode on
dcm = datacursormode(gcf);
set(dcm,'UpdateFcn',@myupdatefcn)

%% --------- Genral Title--------------------------------------------
sgtitle(['Tag ', num2str(TableGPS.Tag(1)),' observ. ',num2str(TableGPS.UniqueSectionCounter(1)),...
    ' starting at ',char(TableGPS.TimeCont(1))])
end

% Save this fcn somewhere on the path
function txt = myupdatefcn(trash,event)
pos = get(event,'Position');
dts = get(event.Target,'Tag');
txt = {dts,...
    ['X: ',num2str(pos(1))],...
    ['Y: ',num2str(pos(2))]};
end

%% plot circ
function plotCirc(circ,x,y)
figure
plot(x,y,'-k')
hold on
for i=1:length(circ.S)
    if isnan(circ.Area(i))
        plot(x(floor(circ.S(i)):ceil(circ.E(i))),y(floor(circ.S(i)):ceil(circ.E(i))),'.r');
    else
        plot(x(floor(circ.S(i)):ceil(circ.E(i))),y(floor(circ.S(i)):ceil(circ.E(i))),'.g');
    end
end
hold off
axis equal
end