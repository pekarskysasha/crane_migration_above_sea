function [soar_phase,soar_glide]=FindSoarAndGlide(x,y,h,t,magx)

W=10; % Window for differences in heights (seconds)
M1=15;
M2=30;
M3=30;
M4=30; % Minimum gliding time (seconds)
M5=60; % Max time between soar and glide (seconds)
M7=-0.5; % Min glide dz rate (m/s)
M8=10; % Max climb while gliding (m)
M9=10; % Max non-downward steps while gliding (seconds)
M6=0.07; % Min mean climb rate
TO_DRAW=true; % Draw results
TO_DRAW=false; % Draw results

magx=fillmissing(magx,'linear');
magx=magx-min(magx);
magx=magx./max(magx);

climbs=GetClimbs(h,W,M1,M2,M3);
soar_phase=climbs;
soar_phase(:,3:4)=0;

if TO_DRAW
    figure
    subplot(3,1,1)
    plot(x,y)
    hold on
    subplot(3,1,2)
    plot(t,h)
    datetick('x')
    axis tight
    hold on
    for i=1:size(soar_phase,1)
        plot(t(soar_phase(i,1):soar_phase(i,2)),h(soar_phase(i,1):soar_phase(i,2)),'-k','linewidth',1);
    end
    
    subplot(3,1,3)
    plot(t,magx,'-k')
    datetick('x')
    axis tight
end


soar_phase=GetCircles(x,y,soar_phase);
soar_phase=GetNonCircles(x,y,magx,soar_phase,M6);
soar_phase=soar_phase(soar_phase(:,3)>0,:);
if TO_DRAW
    subplot(3,1,2)
    for i=1:size(soar_phase,1)
        plot(t(soar_phase(i,1):soar_phase(i,2)),h(soar_phase(i,1):soar_phase(i,2)),'-k','linewidth',2);
    end
end
soar_glide=GetGlides(h,soar_phase,M4,M5,M7,M8,M9);
soar_glide=soar_glide(soar_glide(:,5)>0,:);
if TO_DRAW
    for i=1:size(soar_glide,1)
        subplot(3,1,1)
        plot(x(soar_glide(i,1):soar_glide(i,2)),y(soar_glide(i,1):soar_glide(i,2)),'-k','linewidth',2);
        plot(x(soar_glide(i,3):soar_glide(i,4)),y(soar_glide(i,3):soar_glide(i,4)),'-r','linewidth',2);
        subplot(3,1,2)
        plot(t(soar_glide(i,1):soar_glide(i,2)),h(soar_glide(i,1):soar_glide(i,2)),'-k','linewidth',3);
        plot(t(soar_glide(i,3):soar_glide(i,4)),h(soar_glide(i,3):soar_glide(i,4)),'-r','linewidth',3);
    end
end
end

%%
function soar_glide=GetGlides(h,soar_phase,M4,M5,M7,M8,M9)
soar_glide=soar_phase;
soar_glide(:,3:5)=0;
%  figure
%     plot(h)
%     hold on
for i=1:size(soar_glide,1)
    se=soar_glide(i,2);
    f=se;
    if (i<size(soar_glide,1))
        q=soar_glide(i+1,1)-1;
    else
        q=length(h);
    end
    g=min(q,f+M4);
    while (g<q) && (max(h(g:min(end,g+10)))-h(g)<M8) && (max(imfilter(double(diff(h(g:min(end,g+10)))>=0),ones(M9,1)))<M9)
        g=g+1;
    end
    
    f=find((h(f:(g-1))>=h((f-1):(g-2))) & (h(f:(g-1))>h((f+1):(g))),1)+f-1;
    while (~isempty(f) && (f<=g-M4) && (f<se+M5))
        if mean(diff(h(f:f+M4))>=0)<0.25
            dz=(h((f+M4):g)-h(f))./(M4:(g-f))';
            if min(dz)<=M7
                if ~isempty(g)
                    dz=(h(g)-h(f:(g-M4)))./((g-f):-1:M4)';
    
                    [~,m]=min(dz);
                    m=min(f+m-1,se+M5);
                    while f<m
                        if ((h(f)-h(m))/(m-f)>M7) || (f==m)
                            f=m;
                            while (f>se+4) && (h(f-5)>=h(f)),f=f-5;end
                            while (f>se) && (h(f-1)>=h(f)),f=f-1;end
                            while (f<se+M5) && (h(f+1)==h(f)),f=f+1;end
                            break;
                        else
                            f=f+1;
                        end
                    end
                    soar_glide(i,3:5)=[f,g,soar_phase(i,3)];
                    break;
                end
            end
        end
        f=find((h((f+1):(g-1))>=h(f:(g-2))) & (h((f+1):(g-1))>h((f+2):(g))),1)+f;
    end
    
end
end

%%
function soar_phase=GetNonCircles(x,y,magx,soar_phase,M6)

for i=1:size(soar_phase,1)
    if soar_phase(i,3)==0
        f=soar_phase(i,1);
        g=soar_phase(i,2);
        
        c=0;
        s1=f;
        while (s1<g) && (magx(s1)>magx(s1+1)),s1=s1+1;end
        s=f;
        while (s<s1) && (magx(s)<magx(s+1)),s=s+1;end
        e=find(abs(magx(s+1:g)-magx(s))>0.5,1)+s;
        while ~isempty(e)
            while (e<g) && (abs(magx(e+1)-magx(s))>=abs(magx(e)-magx(s))),e=e+1;end
            %disp(['e-s ',num2str([e-s,s,e])])
            if (e-s<60),c=c+1;end
            s=e;
            e=find(abs(magx(s+1:g)-magx(s))>0.6,1)+s;
        end
        c=floor(c/2);
        if c>=2
            soar_phase(i,3)=2;
            soar_phase(i,4)=c;
        end
        
        %         az=atan2(diff(y(f:g)),diff(x(f:g)));
        %         daz=mod(diff(az),2*pi);
        %         daz(daz>pi)=daz(daz>pi)-2*pi;
        %         nmx=[];
        %         for j=1:length(daz)-10
        %             nmx=[nmx,abs(nanmean(daz(j:j+10)))];
        %         end
        %         %         tmp=[tmp nanmean(nmx)];
        %         %         disp(nanmean(nmx))
        %         [c,nanmean(nmx)>M6]
        %         if nanmean(nmx)>M6
        %             soar_phase(i,3)=2;
        %         end
    end
end


% title(num2str(tmp))

end

%%
function soar_phase=GetCircles(x,y,soar_phase)
for i=1:size(soar_phase,1)
    f=soar_phase(i,1);
    g=soar_phase(i,2);
    pcross=zeros(g-f-14,2);
    for j=f:g-15
        f1=max(1,j-1);
        f2=min(j+30,g);
        pcross(j-f+1,:)=FindCrossing(x(f1:f2),y(f1:f2))+f1;
    end
    pcross=pcross(~isnan(pcross(:,1)),:);
    c=size(pcross,1);
    %disp([i,c])
    if c>=2
        soar_phase(i,3)=1;
        soar_phase(i,4)=c;
    end
end
end

%%
function climbs=GetClimbs(h,W,M1,M2,M3)
dh=h(W+1:end)-h(1:end-W);
f=find(dh>0,1);
climbs=nan(0,2);
while ~isempty(f)
    g=find(dh(f+1:end)<0,1)+f;
    if isempty(g),g=length(dh);end
    if g-f>=M1
        climbs=[climbs;f,g];
    end
    
    f=find(dh(g+1:end)>0,1)+g;
end

i=1;
while i<size(climbs,1)
    if climbs(i+1,1)-climbs(i,2)<=M2
        climbs(i,2)=climbs(i+1,2);
        climbs=climbs([1:i,i+2:end],:);
    else
        i=i+1;
    end
end

for i=1:size(climbs,1)
    m=find(h(climbs(i,1):climbs(i,2))==min(h(climbs(i,1):climbs(i,2))),1,'last');
    climbs(i,1)=climbs(i,1)+m-1;
    [~,m]=max(h(climbs(i,1):(climbs(i,2)+W)));
    climbs(i,2)=climbs(i,1)+m-1;
end
climbs=climbs(climbs(:,2)-climbs(:,1)>=M3,:);


end

function pcross=FindCrossing(x,y)
if ((x(1)==x(2)) && (y(1)==y(2))), pcross=[nan,nan];return;end
while length(x)>2 &&((x(3)==x(2)) && (y(3)==y(2)))
    x=x([1,2,4:end]);
    y=y([1,2,4:end]);
end
xi=x(1:2);
yi=y(1:2);
xj=x(3:end);
yj=y(3:end);
[px,py,pind]=polyxpoly(xi,yi,xj,yj);
if (isempty(px)), pcross=[nan,nan];return;end
while (length(px)>1)
    xj=xj(1:end-1);
    yj=yj(1:end-1);
    [px,py,pind]=polyxpoly(xi,yi,xj,yj);
end
if (isempty(px)), pcross=[nan,nan];return;end
if (xi(1)==xi(2)),ap=(py-yi(1))/diff(yi);else,ap=(px-xi(1))/diff(xi);end
if length(xj)<2
    pcross=[nan,nan];return;
end
bx=(px-xj(pind(2)))./(xj(pind(2)+1)-xj(pind(2)));
by=(py-yj(pind(2)))./(yj(pind(2)+1)-yj(pind(2)));
if (isinf(bx)), bx=nan;end
if (isinf(by)), by=nan;end
if (size(pind,1)~=1), error('k'); end
pcross=[ap-1,pind(2)+nanmean([bx,by])];
end