% main function 
dates = ["0216", "0217", "0218", "0221", "0222", "0223", "0225", "0228", "0301"];
trialDataAll = struct.empty;
for d = 1:length(dates)
    date = dates(d);
    taskfile = dir("../data/preprocessed/2022"+date+ "/*_BEH_*.mat");
    load(taskfile.folder +"/"+ taskfile.name);
    % for i = 1:length(trialData)
    %     trialData(i).reactionTime = behavData(i).t4_reactionTime - behavData(i).t3_goCueTime;
    %     trialData(i).peakSpeed = behavData(i).peakSpeed;
    %     trialData(i).delayTime = behavData(i).t3_goCueTime - behavData(i).t2_targetOnsetTime;
    %     trialData(i).trialNew = trialData(i).trial + d*10000;
    % end
    reactionTimes = [behavData.t4_reactionTime] - [behavData.t3_goCueTime];
    peakSpeeds = [behavData.peakSpeed];

    % reactionTimes and peakSpeeds contain NaNs
    % continue to use the data without NaNs removed
    zscoreRT = nan(1, length(reactionTimes));
    curInds = ~isnan(reactionTimes) & reactionTimes > 250 & reactionTimes < 400;
    X = reactionTimes(curInds);
    % Calculate median and MAD
    X_median = median(X);
    X_mad = mad(X, 1);  % '1' makes it consistent with standard deviation for normally distributed data

    % Calculate Modified Z-score
    X_modZscored = 0.6745 * (X - X_median) / X_mad;  % 0.6745 is approximately the Z-score at the 0.75 quantile of a standard normal distribution
    zscoreRT(curInds) = X_modZscored;

    zscorePS = nan(1, length(peakSpeeds));
    curInds = find(~isnan(peakSpeeds));
    X = peakSpeeds(curInds);
    X_median = median(X);
    X_mad = mad(X, 1);
    X_modZscored = 0.6745 * (X - X_median) / X_mad;
    zscorePS(curInds) = X_modZscored;
    for i = 1:length(behavData)
        behavData(i).zReactionTime = zscoreRT(i);
        behavData(i).zPeakSpeed = zscorePS(i);
        behavData(i).reactionTime = reactionTimes(i);
        behavData(i).delayTime = behavData(i).t3_goCueTime - behavData(i).t2_targetOnsetTime;
    end

    save(taskfile.folder +"/"+ taskfile.name, 'behavData');
end