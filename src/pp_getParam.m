function [dates, rootFolder, axisOutputFileName, processedFileNames, analysisBin, params] = pp_getParam(options)
    arguments
        options.Dataset {mustBeMember(options.Dataset, ["all"])}
        options.timeperiod {mustBeMember(options.timeperiod, ["GC", "HT", "TO", "MO", "SU", "DF", "RF"])}
        options.neuralDataType {mustBeMember(options.neuralDataType, ["", "partial-stitched", "whole-stitched"])}
    end
    dataset = options.Dataset; timeperiod = options.timeperiod; neuralDataType = options.neuralDataType;
    addpath(genpath("./utils/plotting/"));
    if dataset == "all"
        dates = ["0216", "0217", "0218","0221","0222","0223","0225","0228", "0301"];
    end

    if timeperiod == "GC"
        tpName = "GoCue";
    elseif timeperiod == "HT"
        tpName = "TargetHold";
    elseif timeperiod == "TO"
        tpName = "TargetOnset";
    elseif timeperiod == "MO"
        tpName = "MovementOnset";
    elseif timeperiod == "SU"
        tpName = "AfterSuccess";
    elseif timeperiod == "DF"
        tpName = "DelayFailure";
    elseif timeperiod == "RF"
        tpName = "ReachFailure";
    end

    
    rootFolder = "../results/neuralData_" + tpName;
    axisOutputFileName = "../interim/"+neuralDataType+"_"+timeperiod+"_rewardAxis_"+dataset+".mat";
    processedFileNames = [];
    if timeperiod == "GC"
        processedFileName = "_GC_200_600_3456.mat";
    elseif timeperiod == "HT"
        processedFileName = "_HT_200_375_3456.mat";
    elseif timeperiod == "TO"
        processedFileName = "_TO_200_600_3.mat";
    elseif timeperiod == "MO"
        processedFileName = "_GC_200_800_3456.mat";
    elseif timeperiod == "SU"
        processedFileName = "_SU_200_800_34567.mat";
    elseif timeperiod == "DF"
        processedFileName = "_DF_200_2000_311.mat";
    elseif timeperiod == "RF"
        processedFileName = "_RF_200_1000_3412.mat";
    end
    for i = 1:length(dates)
        processedFileNames = [processedFileNames; "../data/processed/"+neuralDataType+"/2022"+dates(i)+processedFileName];
    end

    if timeperiod == "GC"
        analysisBin = (50:250);
    elseif timeperiod == "HT"
        analysisBin = (300:500);
    elseif timeperiod == "TO"
        analysisBin = (400:600);
    elseif timeperiod == "MO"
        analysisBin = nan;
    elseif timeperiod == "SU"
        analysisBin = (200:400);
    elseif timeperiod == "DF"
        analysisBin = (700:900);
    elseif timeperiod == "RF"
        analysisBin = (500:700);
    end

    if timeperiod == "GC"
        Xlim = [-200 600];
    elseif timeperiod == "HT"
        Xlim = [-200 375];
    elseif timeperiod == "TO"
        Xlim = [-200 600];
    elseif timeperiod == "MO"
        Xlim = [-200 800];
    elseif timeperiod == "SU"
        Xlim = [-200 800];
    elseif timeperiod == "DF"
        Xlim = [-200 2000];
    elseif timeperiod == "RF"
        Xlim = [-200 1000];
    end
    params.Xlim = Xlim;

    %% data around GoCue -200~+600 ms
    if timeperiod == "TO"
        state = 3;
    elseif timeperiod == "GC"
        state = 4;
    elseif timeperiod == "MO"
        state = 4;
    elseif timeperiod == "HT"
        state = 6;
    elseif timeperiod == "SU"
        state = 7;
    elseif timeperiod == "DF"
        state = 11;
    elseif timeperiod == "RF"
        state = 12;
    end
    params.state = state;
    params.timeperiod = timeperiod;
    params.dataset = dataset;
    params.tpName = tpName;