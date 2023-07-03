function VS_delay_time(data, rewardLabels, difficultyLabels, delayTimes, label, outputFolder)
    nrewards = 3; rewards = [1 2 3];
    ndifficulties = 2; difficulties = [0 1];
    rewColors = [1 0 0; 1 0.6470 0; 0 0 1]; diffColors = [0 0.447 0.741; 0.466 0.674 0.188];
    direColors = {[1 .5 .5],[.75 .75 .5],[.5 1 .5],[.25 .75 .5],[0 .5 .5],[0.25 0.25 .5],[0.5 0 .5],[0.75 0.25 .5]};
    DiffStyle = ["-", ":"];
    DelayTimes = 300:50:1050; nDelayTimes = length(DelayTimes);

    figure; hold on;

    for j=1:ndifficulties
        for i=1:nrewards
            Y = zeros(nDelayTimes, 1); Yerr = zeros(nDelayTimes, 1);
            for k=1:nDelayTimes
                curInds = rewardLabels == rewards(i) & difficultyLabels == difficulties(j) & delayTimes == DelayTimes(k) & ~isnan(data);
                Y(k) = mean(data(curInds));
                Yerr(k) = std(data(curInds)) / sqrt(sum(curInds));
            end
            h(j) = errorbar(1:nDelayTimes, Y,Yerr, 'Color', rewColors(i, :), 'LineWidth', 2.5,  'LineStyle', DiffStyle(j));
        end
    end
    set(gca, 'fontsize', 20, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
    xlim([0.7 nDelayTimes+.3]); xticks(1:3:nDelayTimes); xticklabels(DelayTimes(1:3:end));
    xlabel("Delay Length"); ylabel(label); legend(h, ["Tiny", "Huge"], Location="best");
    set(gcf,'position',[0,0,550,550]);
    saveas(gcf, outputFolder+"-vs-DelayLength.jpg");
    close all;
end