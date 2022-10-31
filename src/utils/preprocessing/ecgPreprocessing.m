function [ denoisedECG] = ecgPreprocessing(ecg,fs)

denoisedECG = double(ecg);
for s=(1:8) % notch filtering
    d = designfilt('bandstopiir', 'filterOrder', 2, ...
                    'HalfPowerFrequency1', 60*s-1, 'HalfPowerFrequency2', 60*s+1, ...
                    'DesignMethod', 'butter', 'SampleRate', fs);
    denoisedECG = filtfilt(d, denoisedECG);
end
denoisedECG = downsample(denoisedECG, 10);