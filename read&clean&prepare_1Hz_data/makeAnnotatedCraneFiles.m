%% annotate all 1Hz data
if 1
TableGPSAllNew.datetime = dateshift(TableGPSAllNew.TimeCont,'start','hour','nearest');
% Take a block of equal hour (the datetime column), and take the median
% point of that block. Or, take the point closest to the rounded hour.
sections = unique(TableGPSAllNew.UniqueSectionCounter);
annotatedTable = [];
dataWithQuartile = table;
for i = 1:length(sections)
    display(['Now in section ',mat2str(sections(i))]);
    data = TableGPSAllNew(TableGPSAllNew.UniqueSectionCounter==sections(i),:);
    data.Properties.VariableNames{20} = 'lat'; % hard-coded to change column "Interpolated_Lat" to "lat"
    data.Properties.VariableNames{21} = 'lon'; % hard-coded to change column "Interpolated_Lon" to "lon"
    [ds, dsp] = loadECMWFdata(data.datetime(1));
    quartile = data.TimeCont;
    quartile.Minute = round(minute(data.TimeCont)/15)*15;
    quartile.Second = 0;
    data.quartile = quartile;
    qHourBlocks = unique(data.quartile);
    inputTable = table;
    for h = 1:size(qHourBlocks,1)
        [~,timeIndex] = min(abs(data.TimeCont-qHourBlocks(h)));
        inputTable(h,:) = annotateTrack(data(timeIndex,:),ds,dsp);
    end
    annotatedTable = [annotatedTable; inputTable];
    dataWithQuartile = [dataWithQuartile; data];
end
% this loop adds the Individual column to TableGPSAllNew, from StatsTCleaned.
for i = 1:size(StatsTCleaned,1)
    TableGPSAllNew.Individual(TableGPSAllNew.UniqueSectionCounter==StatsTCleaned.UniqueSectionCounter(i)) = StatsTCleaned.Individual(i);
end
TableGPSAllNew.quartile = dataWithQuartile.quartile;
end
TableGPSAllNew_Original = TableGPSAllNew;
tags = unique(TableGPSAllNew_Original.Individual);
for tag = 1:length(tags)
    sections = unique(TableGPSAllNew_Original.UniqueSectionCounter(TableGPSAllNew_Original.Individual==tags(tag)));
    newData = table;
    for i = 1:length(sections)
        display(['Now in section ',mat2str(sections(i))]);
        data = TableGPSAllNew_Original(TableGPSAllNew_Original.UniqueSectionCounter==sections(i),:);
        annotatedData = annotatedTable(annotatedTable.UniqueSectionCounter==sections(i),:);
        for j = 1:size(annotatedData,1)
            temp = data(data.quartile==annotatedData.quartile(j),:);
            for k = 1:size(temp,1)
                newData = [newData; temp(k,:), annotatedData(j,30:52)];
            end
        end
    end
%     eval(['TableGPSAllNew_',mat2str(tags(tag)),' = newData;']);
    TableGPSAllNew = newData;
    save(['1HzGPS_and_ACC_annotated_',mat2str(tags(tag))],'TableGPSAllNew');
end