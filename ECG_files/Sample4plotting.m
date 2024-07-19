function Sample4plotting
% This simple function plots ECGs from collar bone and wrist recording
% during first No-Motion period (Seg1_NoMot). To plot recordings from other
% sessions, simply modify the line below to assign a different field name.
    sampField = 'Seg1_NoMot';
    ecgDir = 'ECGs';
    [datafilenames, L] = FileHandling(ecgDir);
    cleanFlg = true;
    for ii = 1:L
        fileName = datafilenames(ii).name;
        filePath = fullfile(ecgDir, fileName);
        ecgDat = load(filePath);
        % Plotting the entire ECG
%         xx = SamplePlot(ecgDat, ~cleanFlg, sampField);
        % Plotting only clean portions:
        xx = SamplePlot(ecgDat, cleanFlg, sampField, fileName(1:end-4));
        if ~isempty(xx)
            break
        end
    end
end

function xx = SamplePlot(ecgStruct, cleanFlg, sampField, fileName)
    ECG = ecgStruct.ECG_REF.(sampField);
    ECGn = ecgStruct.ECG_Wrist.(sampField);
    segInd = ecgStruct.ecgSegInd.(sampField);
    fs =  ecgStruct.fs;
    if cleanFlg
        ECG = ecgStruct.ECG_REF.(sampField)(segInd(1):segInd(2));
        ECGn = ecgStruct.ECG_Wrist.(sampField)(segInd(1):segInd(2));
    end
    tVec = (0:length(ECG)-1)/fs;
    ax(1) = subplot(2,1,1);
    plot(tVec, ECG)
    title('From Collar Bone (clean)')
    ylabel('Voltage (V)')
    ax(2) = subplot(2,1,2);
    plot(tVec, ECGn)
    title('From Wrist (noisy)')
    ylabel('Voltage (V)')
    xlabel('Time (s)')
    linkaxes(ax, 'x')
    sgtitle(['ECGs from file: ' fileName])
    drawnow
    xx = input('Hit enter to continue or enter 0 to terminate');
    close
end

function [datafilenames, L] = FileHandling(datafolder2)
    datafilenames = dir(datafolder2);
    L = length(datafilenames);
    m = 1;
    for n = 1:L
        if datafilenames(n).bytes == 0
            out(m) = n;
            m = m + 1;
        end
    end
    datafilenames(out) = [];
    L = length(datafilenames);
end
