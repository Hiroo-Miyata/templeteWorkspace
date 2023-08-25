function VS_delay_time_success_rate(failedLabels, rewardLabels, difficultyLabels, delayTimes, options)

    arguments
        failedLabels
        rewardLabels
        difficultyLabels
        delayTimes
        options.Label
        options.OutputFolder
        options.Ylim = [20 100]
        options.DiffLegendPos = "southwest"
        options.RewLegendPos = "south"
    end

    [rewardNames, rewardLegends, rewColors, diffColors, direColors, DiffStyle, DelayTimes, nDelayTimes] = getExperimentConstants();
    rewards = unique(rewardLabels); nrewards = length(rewards);
    difficulties = unique(difficultyLabels); ndifficulties = length(difficulties);
    

    figure; hold on;

    for j=1:ndifficulties
        for i=1:nrewards
            Y = zeros(nDelayTimes, 1);
            for k=1:nDelayTimes
                curInds = rewardLabels == rewards(i) & difficultyLabels == difficulties(j) & delayTimes == DelayTimes(k);
                curInds2 = curInds & failedLabels;
                Y(k) = 100 * (1-sum(curInds2)/sum(curInds));
            end
            plot(Y, 'Color', rewColors(rewards(i), :), 'LineWidth', 2.5,  'LineStyle', DiffStyle(j));
            lh(j) = plot(0, 0, 'Color', 'k', 'LineWidth', 2.5,  'LineStyle', DiffStyle(j));
            lr(i) = plot(0, 0, 'Color', rewColors(rewards(i), :), 'LineWidth', 2.5);
        end
    end
    set(gca, 'fontsize', 16, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
    xlim([0.7 nDelayTimes+.3]); ylim(options.Ylim); xticks(1:3:nDelayTimes); xticklabels(DelayTimes(1:3:end));
    xlabel("Delay Length"); ylabel(options.Label); set(gcf,'position',[0,0,550,550]);
    legend(lh, ["Tiny", "Huge"], Location=options.DiffLegendPos); ah1=axes('position',get(gca,'position'),'visible','off');
    leg2=legend(ah1,lr,rewardLegends, Location=options.RewLegendPos);
    set(gca, 'fontsize', 20, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
    saveas(gcf, options.OutputFolder+"-vs-DelayLength.jpg");
    close all;
end