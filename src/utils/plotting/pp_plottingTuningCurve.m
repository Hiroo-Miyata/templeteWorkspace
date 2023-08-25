function pp_plottingTuningCurve(factorScores, dayLabels, directionLabels, options)

    arguments
        factorScores (:,:) double
        dayLabels (1,:) double
        directionLabels (1,:) double
        options.OutputFolder (1,:) char
    end

    directions = unique(directionLabels); ndirections = length(directions);
    days = unique(dayLabels); ndays = length(days);
    [ntrials, nLVs] = size(factorScores);
    dayColors = [1 .5 .5;.75 .75 .5;.5 1 .5;.25 .75 .5;0 .5 .5;0.25 0.25 .5;0.5 0 .5;0.75 0.25 .5];

    for LV = 1:nLVs
        figure; hold on;
        for day = 1:ndays
            for direction = 1:ndirections
                idx = (dayLabels == days(day)) & (directionLabels == directions(direction));
                meanTuningCurve(day, direction) = mean(factorScores(idx, LV));
            end
        end
        for i=1:ndays
            plot(meanTuningCurve(i, :), 'Color', dayColors(i, :), 'linewidth', 2); hold on;
        end
        xticks(1:ndirections); xticklabels(["0", "45", "90", "135", "180", "225", "270", "315"]);
        xlabel("Direction (deg)"); ylabel("Factor Score");
        title(sprintf("Tuning Curve for LV %d", LV));
        % legendLabel should be Day 1, Day 2, etc.
        legendLabel = cellfun(@(x) sprintf("Day %d", x), num2cell(1:ndays), 'UniformOutput', false);
        legend(legendLabel, 'Location', 'best');
        set(gca, 'fontsize', 20, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
        set(gcf,'position',[0,0,550,550]);
        saveas(gcf, fullfile(options.OutputFolder, sprintf("TuningCurve_LV%d.png", LV)));
        close(gcf);
    end

end

