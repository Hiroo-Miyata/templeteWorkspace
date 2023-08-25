function [rewardNames, rewardLegends, rewColors, diffColors, direColors, DiffStyle, DelayTimes, nDelayTimes] = getExperimentConstants()

    rewardNames = ["S", "M", "L", "J"];
    rewardLegends = ["Small", "Medium", "Large", "Jackpot"];
    rewColors = [1 0 0; 1 0.6470 0; 0 0 1; 0 0 0];
    diffColors = [0 0.447 0.741; 0.466 0.674 0.188];
    direColors = {[1 .5 .5],[.75 .75 .5],[.5 1 .5],[.25 .75 .5],[0 .5 .5],[0.25 0.25 .5],[0.5 0 .5],[0.75 0.25 .5]};
    DiffStyle = ["-", ":"];
    DelayTimes = 250:50:1000;
    nDelayTimes = length(DelayTimes);
end