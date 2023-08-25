function VS_trial(trialNums, Y, rewardLabels, difficultyLabels, options)

    arguments
        trialNums (:, 1) double
        Y (:, 1) double
        rewardLabels (:, 1) double
        difficultyLabels (:, 1) double
        options.Label (1, :) char
        options.OutputFolder (1, :) char
        options.IsMultipleDays (1, 1) logical = false
        options.DayLabels (:, 1) double
        options.FigWidth (1, 1) double = 550
        options.DiffLegendPos = "southwest"
        options.RewLegendPos = "south"
    end

    [rewardNames, rewardLegends, rewColors, diffColors, direColors, DiffStyle, DelayTimes, nDelayTimes] = getExperimentConstants();
    rewards = unique(rewardLabels); nrewards = length(rewards);
    difficulties = unique(difficultyLabels); ndifficulties = length(difficulties);

    if options.IsMultipleDays
        dayLabels = options.DayLabels;
    else
        dayLabels = ones(length(Y), 1);
    end
    dates = unique(dayLabels); ndates = length(dates);

    figure; hold on;
    for d = 1:ndates
        for i = 1:nrewards
            for j = 1:ndifficulties
                curInds = rewardLabels == rewards(i) & difficultyLabels == difficulties(j) & dayLabels == dates(d) & ~isnan(Y);
                smoothedY = movmean(Y(curInds), 5);
                plot(trialNums(curInds), smoothedY, 'Color', rewColors(rewards(i), :), 'LineWidth', 2.5,  'LineStyle', DiffStyle(j));
                lh(j) = plot(0, 0, 'Color', 'k', 'LineWidth', 2.5,  'LineStyle', DiffStyle(j));
                lr(i) = plot(0, 0, 'Color', rewColors(rewards(i), :), 'LineWidth', 2.5);
            end
        end
        % add xline at the end of each day
        xline(trialNums(find(dayLabels == dates(d), 1, 'last')), 'Color', 'k', 'LineWidth', 2.5);
    end

    set(gca, 'fontsize', 20, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
    xlabel("Trial Number"); ylabel(options.Label);
    set(gcf,'position',[0,0,options.FigWidth,550]);
    legend(lh, ["Tiny", "Huge"], Location=options.DiffLegendPos); ah1=axes('position',get(gca,'position'),'visible','off');
    leg2=legend(ah1,lr,rewardLegends, Location=options.RewLegendPos);
    set(gca, 'fontsize', 20, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
    saveas(gcf, options.OutputFolder+"-vs-trial.jpg");
    close all;
end