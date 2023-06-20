clear all
close all
%% Make list of Indeviduals
directory='../DataArchive';
FilesIn=struct2table(dir(directory));
IndevidualsExistTemp=FilesIn.name(FilesIn.bytes>0);

DateStart=datenum({'01/10/18';'01/10/19';'01/10/20'},'dd/mm/yy');
DateEnd=datenum({'01/04/19';'01/04/20';'01/04/21'},'dd/mm/yy');
seasomChoose=[{'Fall','Spring'}];
load CraneColors
CraneColors=[CraneColors; CraneColors; CraneColors];

% create variables
Tag=[];
Individual=[];
TimestartT=[];
TimeEndT=[];
TimeOfBurstT=[];
BatteryDranageT=[];
SpeedVsBat=[];
BatEnd=[];
BatStart=[];
MeadianInterval=[];
Season=[];
Year=[];

TagArea=[];
IndividualArea=[];
TimestartArea=[];
TimeEndArea=[];
TAreaAll=[];

figure
title(['fall 2018']);
hold on
figure
title(['spring 2019']);
hold on
figure
title(['fall 2019']);
hold on
figure
title(['spring 2020']);
hold on
figure
title(['fall 2020']);
hold on
figure
title(['spring 2021']);
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-if you want to add to older version load
SratsT=[];
load 1HertzStats_all
add=0;
if ~isempty(SratsT)
    add=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for kk=1:length(IndevidualsExistTemp)
    load([directory,'\',char(IndevidualsExistTemp(kk))]);
    % indevidual
    dd=char(IndevidualsExistTemp(kk));
    indev=str2double(dd(1:end-4));
    Secnd=0;
    for ii=1:length(DateStart)
        Relevantdata=floor(TagMAT(:,3))>=DateStart(ii) & floor(TagMAT(:,3))<=DateEnd(ii) &...
            ~isnan(TagMAT(:,8))  & TagMAT(:,5)~=0;
        if sum(Relevantdata)==0 | indev<5000 % no e-Obs tags
            continue
        end
        % check intervals
        mldate=TagMAT(Relevantdata,3);
        Interv=diff(mldate)*86400; % convert intervals to sec
        % find 1 Hz
        HZ1=Interv<=7 & Interv>0;
        if sum(HZ1)==0 % if no 1 HZ
            continue
        end
        % find start and end of burst
        INDj=diff(HZ1);
        
        if INDj(1)==-1
            INDj=INDj(2:end);
        end
        startBurst=find(INDj==1)+1;
        endBurst=find(INDj==-1);
        if endBurst(end)<startBurst(end)
            startBurst=startBurst(1:end-1);
        end
        IndexBurstTemp=[startBurst,endBurst];
        % find mistake stops
        MsInt=[nan; IndexBurstTemp(2:end,1)-IndexBurstTemp(1:end-1,2)];
        IndexBurst=[IndexBurstTemp(1,1),IndexBurstTemp(1,2)];
        i=2;
        while i<length(MsInt)
            if MsInt(i)<3  | isnan(MsInt(i)) % if the interval small close the gap
                IndexBurst(end,2)=IndexBurstTemp(i,2);
                i=i+1;
            else
                IndexBurst=[IndexBurst; [IndexBurstTemp(i,1),IndexBurstTemp(i,2)]];
                i=i+1;
            end
        end
        
        % eraze very short bursts (less then 10 min)
        I=find([IndexBurst(:,2)-IndexBurst(:,1)]<600);
        IndexBurst(I,:)=[];
        if isempty(IndexBurst)
            continue
        end
        lon=TagMAT(Relevantdata,4);
        lat=TagMAT(Relevantdata,5);
        speed=TagMAT(Relevantdata,8);
        ele=TagMAT(Relevantdata,7);
        bat=TagMAT(Relevantdata,10);
        tag=TagMAT(1,1);
        
        S=[];
        E=[];
        TArea=[];
        seas=[];
        yy=[];
        % Check what happened that the burst ended
        for n=1:length(IndexBurst(:,1))
            T=(mldate(IndexBurst(n,2))-mldate(IndexBurst(n,1)))*86400/60; %time in Minutes
            Timestart(n,1)=mldate(IndexBurst(n,1));
            TimeEnd(n,1)=mldate(IndexBurst(n,2));
            TimeOfBurst(n,1)=T;
            B=(bat(IndexBurst(n,2))-bat(IndexBurst(n,1)));
            BatteryDranage(n,1)=abs(B)/T;
            BatE(n,1)=bat(IndexBurst(n,2));
            BatS(n,1)=bat(IndexBurst(n,1));
            IntervalM(n,1)=median(Interv(IndexBurst(n,1):(IndexBurst(n,2))));
            
            if speed(IndexBurst(n,2))<2.8 & bat(IndexBurst(n,2))>=50
                SpVsBat(n,1)=1; % speed is the reason to stop
            elseif speed(IndexBurst(n,2))>=2.8 & bat(IndexBurst(n,2))<50
                SpVsBat(n,1)=0; % bat is the reason to stop
            elseif speed(IndexBurst(n,2))<2.8 & bat(IndexBurst(n,2))<50
                SpVsBat(n,1)=2; % both
            elseif speed(IndexBurst(n,2))>=2.8 & bat(IndexBurst(n,2))>=50
                SpVsBat(n,1)=3; % Geofence the reason to stop
            end
            %% find location
            lonBurst=lon(IndexBurst(n,1):IndexBurst(n,2));
            latBurst=lat(IndexBurst(n,1):IndexBurst(n,2));
            mldateBurst=mldate(IndexBurst(n,1):IndexBurst(n,2));
            
            LonB=[34.5,35.7];
            LatB=[34.4,36.5];
            InArea=lonBurst>LonB(1) & lonBurst<LonB(2) & latBurst>LatB(1) & latBurst<LatB(2);
            if sum(InArea)>100
                S=[S; min(mldateBurst(InArea))];
                E=[E; max(mldateBurst(InArea))];
                TArea=[TArea; (max(mldateBurst(InArea))-min(mldateBurst(InArea)))*86400/60]; %time in Minutes
            end
            
            % year and season
            if month(mldate(IndexBurst(n,1)))>5 & month(mldate(IndexBurst(n,1)))<12 % fall
                season=1;
            elseif month(mldate(IndexBurst(n,1)))>=1 % spring
                season=2;
            end
            y=year(mldate(IndexBurst(n,1)));
            seas=[seas; seasomChoose(season)];
            yy=[yy; y];
            %% plot only if longer than one hour
            %call the relavant plot
            Legend=num2str(indev);
            if T>60 & median(Interv(IndexBurst(n,1):(IndexBurst(n,2))))<4
                % choose what figure to call
                if  season==1 & y==2018
                    figure(1)
                elseif season==2 & y==2019
                    figure(2)
                elseif season==1 & y==2019
                    figure(3)
                elseif season==2 & y==2020
                    figure(4)
                elseif season==1 & y==2020
                    figure(5)
                elseif season==2 & y==2021
                    figure(6)
                end
                if Secnd==1
                    plot(lon(IndexBurst(n,1):IndexBurst(n,2)),lat(IndexBurst(n,1):IndexBurst(n,2)),'-',...
                        'LineWidth',2,'Color',CraneColors(kk,:),'HandleVisibility','off') %line
                else
                    plot(lon(IndexBurst(n,1):IndexBurst(n,2)),lat(IndexBurst(n,1):IndexBurst(n,2)),'-',...
                        'LineWidth',2,'Color',CraneColors(kk,:),'DisplayName', Legend) %line
                    Secnd=1;
                end
                hold on
            end
            
        end
        %% if we are adding to excsting one then skip if ut akready excists
        index=[1:1:length(Timestart)]';
        if add==1
            RemoveInd=[];
            for t=1:length(Timestart)
                in=SratsT.Individual==indev & datenum(SratsT.TimestartV)==Timestart(t);
                if sum(in)>0
                    RemoveInd=[RemoveInd; t];
                end
            end
        end
        index(RemoveInd)=[];
        % if bo aading is needed
        if isempty(index)
            clear  Timestart TimeEnd TimeOfBurst BatteryDranage SpVsBat BatE BatS IntervalM
            continue;
        end
        %% For area
        IndividualArea=[TagArea; ones(length(S),1)*indev];
        TagArea=[TagArea; ones(length(S),1)*tag];
        TimestartArea=[TimestartArea; S];
        TimeEndArea=[TimeEndArea; E];
        TAreaAll=[TAreaAll; TArea];
        
        Tag=[Tag; ones(length(Timestart(index)),1)*tag];
        Individual=[Individual; ones(length(Timestart(index)),1)*indev];
        TimestartT=[TimestartT; Timestart(index)];
        TimeEndT=[TimeEndT; TimeEnd(index)];
        TimeOfBurstT=[TimeOfBurstT; TimeOfBurst(index)];
        BatteryDranageT=[BatteryDranageT; BatteryDranage(index)];
        SpeedVsBat=[SpeedVsBat; SpVsBat(index)];
        BatEnd=[BatEnd; BatE(index)];
        BatStart=[BatStart; BatS(index)];
        MeadianInterval=[MeadianInterval; IntervalM(index)];
        Season=[Season; seas(index)];
        Year=[Year; yy(index)];
        
        
        clear  Timestart TimeEnd TimeOfBurst BatteryDranage SpVsBat BatE BatS IntervalM
        
    end
end
for f=1:6
    figure(f)
    %legend show
    plot_google_map('MapType','satellite')
    plot([LonB(1),LonB(1)],[LatB(1),LatB(2)],'-r')
    plot([LonB(2),LonB(2)],[LatB(1),LatB(2)],'-r')
    plot([LonB(1),LonB(2)],[LatB(1),LatB(1)],'-r')
    plot([LonB(1),LonB(2)],[LatB(2),LatB(2)],'-r')
end
TimestartVArea=datetime(TimestartArea,'ConvertFrom','datenum','format','y-MM-dd HH:mm:ss');
TimeEndVArea=datetime(TimeEndArea,'ConvertFrom','datenum','format','y-MM-dd HH:mm:ss');
SratsArea=table(IndividualArea,TagArea,TimestartVArea,TimeEndVArea,TAreaAll,....
    'VariableNames',{'Individual','Tag','TimestartV','TimeEndV','TimeOfBurstT'});
save('1HertzStatsMediterranean','SratsArea')
% save the table before adding the new data
SratsTOld=SratsT;
% add new data
TimestartV=datetime(TimestartT,'ConvertFrom','datenum','format','y-MM-dd HH:mm:ss');
TimeEndV=datetime(TimeEndT,'ConvertFrom','datenum','format','y-MM-dd HH:mm:ss');
if add==1
    temp=table(Individual,Tag,Year,Season,...
        TimestartV,TimeEndV,TimeOfBurstT,MeadianInterval,BatteryDranageT,SpeedVsBat,BatStart,BatEnd);
    SratsT=[SratsT; temp];
else
    SratsT=table(Individual,Tag,Year,Season,...
        TimestartV,TimeEndV,TimeOfBurstT,MeadianInterval,BatteryDranageT,SpeedVsBat,BatStart,BatEnd);
end
save('1HertzStats_all_summer2021','SratsT')

%---- for downloading files from ornirela website (SAVE BEFORE) ---------------------------------------
if add==1
    SratsT=SratsT(height(SratsTOld)+1:end,:);
end
SratsT=SratsT(SratsT.MeadianInterval<4,:);
T=unique(SratsT.Tag);
Down=[];
for i=1:length(T)
    st=floor(min(datenum(SratsT.TimestartV(SratsT.Tag==T(i)))));
    ed=floor(max(datenum(SratsT.TimeEndV(SratsT.Tag==T(i)))));
    Down=[Down; [T(i),st,ed]];
end
TimestartVd=datetime(Down(:,2),'ConvertFrom','datenum','format','y-MM-dd')-1;
TimeEndVd=datetime(Down(:,3),'ConvertFrom','datenum','format','y-MM-dd')+1;
DownT=table(Down(:,1),TimestartVd,TimeEndVd);
