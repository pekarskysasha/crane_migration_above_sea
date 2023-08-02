function [wang_calc,wvel_calc]=WindCalcultedThermals(Indexes,TableGPSCont)
% calculations from: Using High-Resolution GPS Tracking Data of Bird Flight for Meteorological Observations
% Treep et al., 2016, Bulletin of the American Meteorological Society
for i=1:size(Indexes,1)
    % calculate time in thermal (seconds)
    timeInThermal=Indexes(i,2)-Indexes(i,1);
    if timeInThermal>72 % calculate only for thermals of at least 72 s 
        lat=TableGPSCont.Interpolated_Lat(Indexes(i,1):Indexes(i,2));
        lon=TableGPSCont.Interpolated_Lon(Indexes(i,1):Indexes(i,2));
        % linear regression through all GPS points
        mdl = fitlm(lon,lat);
        Cof=mdl.Coefficients;
        % calculate the straight line
        X=lon;
        Y=Cof.Estimate(1)+Cof.Estimate(2)*X;
        % net horizontal displacement
        NetDist=distance(Y(end),X(end),Y(1),X(1),wgs84Ellipsoid);
        % The net horizontal displacement of the bird divided by the time span of each circling bout yields an estimate for wind speed
        wvel_calc(i,1)=NetDist/timeInThermal;
        wang_calc(i,1)=atan2(Y(end)-Y(1),X(end)-X(1));
        winddir_cal = atan2d(Y(end)-Y(1),X(end)-X(1)) + 360*(Y(end)-Y(1)<0);
    else
        wang_calc(i,1)=nan;
        wvel_calc(i,1)=nan;
    end
end
end