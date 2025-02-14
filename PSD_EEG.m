dataPaths = {
    'C:\Users\imanf\Downloads\WOS1\chb01_02.edf',
    'C:\Users\imanf\Downloads\WS1\chb01_03.edf',
    'C:\Users\imanf\Downloads\WS1\chb01_04.edf',
    'C:\Users\imanf\Downloads\WS1\chb01_15.edf',
    'C:\Users\imanf\Downloads\WS1\chb01_18.edf',
    'C:\Users\imanf\Downloads\WS1\chb01_16.edf',
    'C:\Users\imanf\Downloads\WS1\chb01_26.edf'
};

EEGData = cell(1, length(dataPaths));
for i = 1:length(dataPaths)
    EEGData{i} = edfread(dataPaths{i});
end
numericData = cell(1, length(EEGData));
for i = 1:length(EEGData)
    if istimetable(EEGData{i})
        numericData{i} = EEGData{i}.Variables;
    else
        numericData{i} = EEGData{i};
    end
end
f = cell(1, length(EEGData)); 
PSD = cell(1, length(EEGData));
for i = 1:length(EEGData)

    [pxx, freq] = pwelch(numericData{i}, [], [], [], 256); 
    PSD{i} = 10*log10(pxx); 
    f{i} = freq; 
disp(['Signal ' num2str(i) ' has unique PSD values: ' num2str(length(unique(PSD{i})))]);

end

for i = 1:length(PSD)
    figure;
    plot(f{i}, PSD{i});
    grid on; 
    title(sprintf('Power Spectral Density of Signal %d', i));
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
end
