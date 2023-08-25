function VS_time(data, rewardLabels, difficultyLabels, options)
    
    arguments
        data (:,:) double
        rewardLabels
        difficultyLabels
        options.Label
        options.OutputFolder
        options.Xlim (1,2) double = nan
        options.Ylim = nan
        options.Timeperiod = nan
        options.DiffLegendPos = "southwest"
        options.RewLegendPos = "south"
    end

    [rewardNames, rewardLegends, rewColors, diffColors, direColors, DiffStyle, DelayTimes, nDelayTimes] = getExperimentConstants();
    rewards = unique(rewardLabels); nrewards = length(rewards);
    difficulties = unique(difficultyLabels); ndifficulties = length(difficulties);

    [timeBin, ntrials] = size(data);
    if isnan(options.Xlim)
        options.Xlim = [0, timeBin];
        Xticks = [0, timeBin];
    else
        if options.Xlim(1) > 0
            Xticks = [options.Xlim(1), options.Xlim(2)];
            Xticklabels = [num2str(options.Xlim(1)), num2str(options.Xlim(2))];
        else
            Xticks = [options.Xlim(1), 0, options.Xlim(2)];
            Xticklabels = [num2str(options.Xlim(1)), options.Timeperiod, num2str(options.Xlim(2))];
        end
    end


    figure; hold on;
    for j=1:ndifficulties
        for i=1:nrewards
            curInds = rewardLabels == rewards(i) & difficultyLabels == difficulties(j);
            Y = mean(data(:, curInds), 2)';
            % calculate 95%CI
            Yerr = 1.96*std(data(:, curInds), 0, 2)/sqrt(sum(curInds)); Yerr = Yerr';

            % plot
            plot(options.Xlim(1):options.Xlim(2), Y, 'Color', rewColors(rewards(i), :), 'LineWidth', 2.5,  'LineStyle', DiffStyle(j));
            lh(j) = plot(0, 0, 'Color', 'k', 'LineWidth', 2.5,  'LineStyle', DiffStyle(j));
            lr(i) = plot(0, 0, 'Color', rewColors(rewards(i), :), 'LineWidth', 2.5);
            fill([options.Xlim(1):options.Xlim(2), options.Xlim(2):-1:options.Xlim(1)], [Y+Yerr, Y(end:-1:1)-Yerr(end:-1:1)], rewColors(rewards(i), :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end
    end
    set(gca, 'fontsize', 20, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
    xlim(options.Xlim); xticks(Xticks); xticklabels(Xticklabels);
    ylabel(options.Label); ylim(options.Ylim);
    set(gcf,'position',[0,0,550,550]);
    legend(lh, ["Tiny", "Huge"], Location=options.DiffLegendPos); ah1=axes('position',get(gca,'position'),'visible','off');
    leg2=legend(ah1,lr,rewardLegends, Location=options.RewLegendPos);
    set(gca, 'fontsize', 20, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
    saveas(gcf, options.OutputFolder+"-vs-time.jpg"); close all;
end