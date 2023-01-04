function plotHeatMap(firingRateAll, badTrials, badUnits, outputFolder)
    
% this function plots the heatmap of firing rates of every good and bad neuron in every good and bad trial
% and help to make sure the p1_spikePreprocessing.m works well.


%% pre plotting: check the size of input
if size(firingRateAll, 1) ~= length(badUnits)
    error("INPUT ERROR: the size of neuralData.channel and batUnits are expected same; \n")
end
if size(firingRateAll, 2) ~= length(badTrials)
    error("INPUT ERROR: the length of neuralData and batTrials are expected same; \n")
end

%% plot heatmap
% plot 3 condition bad trials, bad units, good trials and good units
figure;colormap("turbo");imagesc(firingRateAll(:, badTrials));colorbar;
title("Heap Map at bad trials"); xlabel("Trials"); ylabel("Neuron Channels");
saveas(gcf, outputFolder+"heatmap-badtrials.jpg");
figure;colormap("turbo");imagesc(firingRateAll(badUnits, ~badTrials));colorbar;
title("Heap Map at bad Units"); xlabel("Trials"); ylabel("Neuron Channels");
saveas(gcf, outputFolder+"heatmap-badunits.jpg");
figure;colormap("turbo");imagesc(firingRateAll(~badUnits, ~badTrials));colorbar;
title("Heap Map at good Units and trials"); xlabel("Trials"); ylabel("Neuron Channels");
saveas(gcf, outputFolder+"heatmap-goodtrial-goodunits.jpg");