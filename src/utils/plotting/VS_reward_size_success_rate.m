function VS_reward_size_success_rate(failedLabels, rewardLabels, difficultyLabels, options)

    arguments
        failedLabels
        rewardLabels
        difficultyLabels
        options.Label
        options.OutputFolder
        options.Ylim = [20 100]
    end


    [rewardNames, rewardLegends, rewColors, diffColors, direColors, DiffStyle, DelayTimes, nDelayTimes] = getExperimentConstants();
    rewards = unique(rewardLabels); nrewards = length(rewards);
    difficulties = unique(difficultyLabels); ndifficulties = length(difficulties);

    figure; hold on;
    for j=1:ndifficulties
        Y = zeros(nrewards,1);
        for i=1:nrewards
            curInds = rewardLabels == rewards(i) & difficultyLabels == difficulties(j);
            curInds2 = curInds & failedLabels;
            Y(i) = 100 * (1-sum(curInds2)/sum(curInds));
        end
        plot(Y, 'Color', diffColors(j, :), 'LineWidth', 2.5);
    end
    set(gca, 'fontsize', 16, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
    ylim(options.Ylim);
    xlim([0.7, nrewards+0.3]); xticks(1:nrewards); xticklabels(rewardNames);
    ylabel(options.Label); legend(["Tiny", "Huge"], Location="best"); set(gcf,'position',[0,0,550,550]);
    saveas(gcf, options.OutputFolder+"-vs-Reward.jpg"); close all;
end