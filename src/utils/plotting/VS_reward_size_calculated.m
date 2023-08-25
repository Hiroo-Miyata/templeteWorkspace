function VS_reward_size_calculated(Y, Yerr, options)

    arguments
        Y (3, 2) double
        Yerr = nan
        options.Label (1, :) char
        options.OutputFolder (1, :) char
    end

    label = options.Label; outputFolder = options.OutputFolder;

    [rewardNames, rewardLegends, rewColors, diffColors, direColors, DiffStyle, DelayTimes, nDelayTimes] = getExperimentConstants();
    rewards = unique([1 2 3]); nrewards = length(rewards);
    difficulties = unique([0 1]); ndifficulties = length(difficulties);

    figure; hold on;
    for j=1:ndifficulties
        if isnan(Yerr)
            plot(1:nrewards, Y(:, j), Color=diffColors(j, :), LineWidth=2);
        else
            errorbar(1:nrewards, Y(:, j),Yerr(:, j), Color=diffColors(j, :), LineWidth=3);
        end
    end
    set(gca, 'fontsize', 20, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
    xlim([0.7 3.3]); xticks(1:3); xticklabels(["S", "M", "L"]); ylabel(label); legend(["Tiny", "Huge"], Location="best");
    set(gcf,'position',[0,0,550,550]);
    saveas(gcf, outputFolder+"-vs-Reward.jpg"); close all;
end