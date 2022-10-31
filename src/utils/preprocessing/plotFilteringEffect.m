function plotFilteringEffect(before, after, fs, mainTitle, outputFolder, extensions)

subplot(1,2,1)
plot((1:length(before))/fs, before); hold on;
plot((1:length(after))/fs, after); hold off;
xlabel('EMG voltage (a.u.)');
ylabel('Time (s)');
legend(["Before", "After"])
title("Time Series EMG Voltage");
subplot(1,2,2)
pwelch(before, 30*fs, [], [], fs); axis([0 500 -inf inf]); hold on;
pwelch(after, 30*fs, [], [], fs); axis([0 500 -inf inf]); hold off;
legend(["Before", "After"])
title("30s Window Welch Power Spectral Density Estimate");
sgtitle(mainTitle)
set(gcf, "position", [10 10 1600 400])

for e=(1:length(extensions))
    saveas(gcf, outputFolder + mainTitle + "." + extensions(e));
    close all
end