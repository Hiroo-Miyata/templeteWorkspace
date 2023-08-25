%%% This script takes synchronized data and preprocesses all data streams.
%%% It saves 8 separate files for a session:
%%% - TSK: task parameter data (e.g., dir/rew labels, state trans., etc.)
%%% - KIN: kinematic data for hand position
%%% - EMG: muscle data
%%% - ECG: heart rate information
%%% - EYE: eye gaze and pupil data
%%% - NER: neural data
%%% - BEH: single trial behavioral metrics (e.g., peak speed, saccade RT)
%%% - taskInfo: metrics relevant across the entire session (i.e., task
%%% parameters, mean waveforms, filter information)
%%%
%%% Aside from taskInfo, all of these are structures with length ntrials
%%% that have fields with values corresponding to each trial. Aside from
%%% taskInfo/BEH, all field values within a trial have length ntime, where
%%% ntime is the number of milliseconds in that trial.
%%%
%%% It's recommended you run each of the preprocessing steps in sections to
%%% evaluate output as you go.
%%%
%%% Adam Smoulder, 11/2/2022
clear all; close all;
dbstop if error

dates = ["20220216", "20220217", "20220218", "20220221", "20220222", "20220223", "20220225", "20220228", "20220301"];

for day = 1:length(dates)
date = dates(day);
dataFolder = "..\..\..\data\";
tmpFileInfo = dir(dataFolder+'synchinterim/*'+date+'*.mat');
load(tmpFileInfo.folder +"/"+ tmpFileInfo.name);

% Edit folder/filenames and file tags (what they start with) here
tmpFileInfo = dir(dataFolder+'synchinterim/*'+date+'*.mat');
spikeAndAnalogDataFile = tmpFileInfo.folder +"/"+ tmpFileInfo.name;
tmpFileInfo = dir(dataFolder+'/raw/TaskData/*'+date+'*');
taskDataFolderList = tmpFileInfo.folder +"/"+ tmpFileInfo.name+'\';
functionFolder = ['..\']; % This is where the functions and classes needed are, at least beneath it
outputFolder = ['..\..\..\data\preprocessed\'+date + '\'] ; % where to save files


clear dataFolder tmpFileInfo

% Things that shouldn't change over days
trialDataTag = 'Trial0'; % what each trial data file starts with
rewardNames = {'Small','Medium','Large','Jackpot'}; % Searches target names for these

% Go to working folder
addpath(functionFolder)

% % Adam uses this to remove paths to other stuff
% rmpath(genpath('D:\AdamMatlab'))
% addpath(genpath(pwd))

%% Load the data
disp('Loading data')
dateString = taskDataFolderList{1}(end-15:end-1); % show it for confirmation
load([taskDataFolderList{1} 'taskInfo.mat'])
load(spikeAndAnalogDataFile)
disp('Loaded synchronized data')

%% Do an initial processing of the behavior and set up the file name to save
[trialData,taskInfo,kinematicData] = pp_preprocessTaskAndKinematicData(taskDataFolderList,rewardNames,synchInfo);
procTimeString = datestr(now,'yyyymm');
saveFilenameBase = ['Rocky' ...
    taskDataFolderList{1}(end-15:end-8) ...
    '_Trials' num2str(trialData(1).trial) ...
    '_' num2str(trialData(end).trial) '_preprocessed_' ...
    ];
save(outputFolder+saveFilenameBase+'TSK_'+procTimeString,'trialData', 'taskInfo')
save(outputFolder+saveFilenameBase+'KIN_'+procTimeString,'kinematicData')
disp('Saved task and kinematic data')

%% Attach spikes to trialData
plotoutputFolder = outputFolder+'neuralPreprocessingPlot/';
mkdir(plotoutputFolder)
[trialData,taskInfo,neuralData] = pp_spikeAlignmentToTrial(taskAlignedSpiketimes,synchInfo,trialData,taskInfo);
[neuralData, spikeInfo] = p1_spikePreprocessing(trialData, neuralData,1,0.35, plotoutputFolder);
save(outputFolder+saveFilenameBase+'NER_'+procTimeString,'neuralData', "spikeInfo")
disp('Saved neural data')

%% Preprocess EMG and ECG data
% TO-DO: put in EMG processing here
% outputFormat = ["jpg"]; %"jpg", "svg", "fig"
% plotoutputFolder = outputFolder+'emgConfirmPlot/';
% mkdir(plotoutputFolder)
% [EMGData, EMGMetrics, ECGData] = pp_emgPreprocess(trialData, analogData, outputFormat, plotoutputFolder);
% save(outputFolder+saveFilenameBase+'EMG_'+procTimeString, 'EMGData', 'EMGMetrics')
% save(outputFolder+saveFilenameBase+'ECG_'+procTimeString, 'ECGData')
% disp('Saved EMG data')

%% Preprocess eye data
% TO-DO: put in Eye data processing here
% outputFormat = ["jpg"]; %"jpg", "svg", "fig"
% plotoutputFolder = outputFolder+'pupilPlot/';
% mkdir(plotoutputFolder)
% pupilData = pp_pupilPreprocess(trialData, analogData, outputFormat, plotoutputFolder);
% save(outputFolder+saveFilenameBase+'EYE_'+procTimeString, 'pupilData')

%% Assess various data streams for behavioral metrics (i.e., peak speed,
% saccade reaction time, average heart rate, etc.)
% TO-DO: put in behavioral processing here
% save([outputFolder saveFilenameBase 'BEH_' procTimeString],'behaviorData')

disp('Completed preprocessing')

end




