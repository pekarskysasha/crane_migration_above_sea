function inputTable = annotateTrack(inputTable, ds, dsp)
    % Usage: inputTable = annotateTrack(inputTable, ds, dsp)
    %
    % annotateTrack annotates inputTable using variables that are in ds and dsp.
    % Currently, variables to annotate are hard-coded into the function.
    % Later we can add the option of taking all the variables directly from
    % ds and dsp.
    % 
    % inputTable: a table of location data, with the following columns:
    % datetime, lat, lon.
    %
    % datetime is a datetime variable that is the ECMWF timestamp nearest to the 
    % real timestamp of this location point (that is, the nearest round hour). 
    % For example, if the real timestamp is '10-Oct-2018 15:12:31', then the 
    % corresponding time is '10-Oct-2018 15:00:00'.
    %
    % lat and lon are the latitude and longitude of the location point.
    if sum(contains(inputTable.Properties.VariableNames,'datetime')) == 0 || sum(contains(inputTable.Properties.VariableNames,'lat')) == 0 || sum(contains(inputTable.Properties.VariableNames,'lon')) == 0
        disp('inputTable not formatted correctly, please use variables "datetime", "lat", "lon".')
        return
    end
    % first create the structures that hold all ECMWF variables of
    % interest, including the lat/lon/time vectors of each nc file
    dsVariables = {'t2m', 'sst', 'blh', 'cbh', 'cape', 'cin', 'msshf', 'tcc', 'msl'};
    for i = 1:length(dsVariables)
        eval(['dsData.',str2mat(dsVariables(i)),'=ds.geovariable(dsVariables(i));']);
    end
    dspVariables = {'u', 'v', 't'};
    dspLevels = dsp.geovariable('level');
    dspLevels = dspLevels(:);
    for i = 1:length(dspVariables)
        eval(['dsData.',str2mat(dspVariables(i)),'=dsp.geovariable(dspVariables(i));']);
    end
    ECvars = fieldnames(dsData);
    for i = 1:size(inputTable,1)
        for j = 1:length(ECvars)
            v = eval(['dsData.',str2mat(ECvars(j))]);
            [~,timeIndex] = min(abs(v.gettimedata(1)-datenum(inputTable.datetime(i))));
            [~,latIndex] = min(abs(v.getlatdata(1)-inputTable.lat(i)));
            [~,lonIndex] = min(abs(v.getlondata(1)-inputTable.lon(i)));
            if isfield(v.getaxesorder,'level') % level axis exists
                if isfield(v.getaxesorder,'expver') % expver axis exists
                    torun = ['d = v(timeIndex,1,level,latIndex,lonIndex);'];
                else
                    torun = ['d = v(timeIndex,level,latIndex,lonIndex);'];
                end
                for level = 1:length(dspLevels)
                    eval(torun);
                    if contains(v.attribute('long_name'),'emperature') && v.attribute('units') == 'K' % if temperature in Kelvin, convert to Celsius
                        d = d - 273.15;
                    end
                    eval(['inputTable.',str2mat(ECvars(j)),int2str(dspLevels(level)),'(i) = d;']);
                end
            else % no level axis
                if isfield(v.getaxesorder,'expver') % expver axis exists
                    d = v(timeIndex,1,latIndex,lonIndex);
                else % only time,lat,lon axes
                    d = v(timeIndex,latIndex,lonIndex);
                end
                if contains(v.attribute('long_name'),'temperature') && v.attribute('units') == 'K' % if temperature in Kelvin, convert to Celsius
                    d = d - 273.15;
                end
                eval(['inputTable.',str2mat(ECvars(j)),'(i) = d;']);
            end
        end
        inputTable.sstDiff(i) = inputTable.sst(i) - inputTable.t2m(i);
        for level = 1:length(dspLevels)
            [windDirection, windSpeed] = eval(['calcWind(inputTable.u',int2str(dspLevels(level)),'(i), inputTable.v',int2str(dspLevels(level)),'(i));']);
            eval(['inputTable.windDirection',int2str(dspLevels(level)),'(i) = windDirection;']);
            eval(['inputTable.windSpeed',int2str(dspLevels(level)),'(i) = windSpeed;']);
            if i < size(inputTable,1) % can't calculate heading if this is the last point
                inputTable.movementHeading(i) = calcWind(inputTable.lon(i+1)-inputTable.lon(i),inputTable.lat(i+1)-inputTable.lat(i));
                eval(['inputTable.windSupport',int2str(dspLevels(level)),'(i) = inputTable.windDirection',int2str(dspLevels(level)),'(i) - inputTable.movementHeading(i);']);
            else
                inputTable.movementHeading(i) = NaN;
                eval(['inputTable.windSupport',int2str(dspLevels(level)),'(i) = NaN;']);
            end
        end
    end
end

function [windDir, windSpeed] = calcWind(u,v)
    windSpeed = sqrt(u.^2 + v.^2);
    windDir = atand(abs(u)/abs(v));
    if u < 0 && v < 0
        windDir = windDir + 180;
    elseif u < 0 && v > 0
        windDir = 360 - windDir;
    elseif u > 0 && v > 0
        windDir = windDir; % this is unnecessary but easier to understand
    elseif u > 0 && v < 0
        windDir = 180 - windDir;
    end
end

