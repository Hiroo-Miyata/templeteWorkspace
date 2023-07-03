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
    end

    nrewards = 3; rewards = [1 2 3];
    ndifficulties = 2; difficulties = [0 1];
    rewColors = [1 0 0; 1 0.6470 0; 0 0 1]; diffColors = [0 0.447 0.741; 0.466 0.674 0.188];
    direColors = {[1 .5 .5],[.75 .75 .5],[.5 1 .5],[.25 .75 .5],[0 .5 .5],[0.25 0.25 .5],[0.5 0 .5],[0.75 0.25 .5]};
    DiffStyle = ["-", ":"];

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
                smoothedY = movmean(Y(curInds), 50);
                h(j) = plot(trialNums(curInds), smoothedY, 'Color', rewColors(i, :), 'LineWidth', 2.5,  'LineStyle', DiffStyle(j));
            end
        end
        % add xline at the end of each day
        xline(trialNums(find(dayLabels == dates(d), 1, 'last')), 'Color', 'k', 'LineWidth', 2.5);
    end

    set(gca, 'fontsize', 20, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
    xlabel("Trial Number"); ylabel(options.Label); legend(h, ["Tiny", "Huge"], Location="best");
    set(gcf,'position',[0,0,options.FigWidth,550]);
    saveas(gcf, options.OutputFolder+"-vs-trial.jpg");
    close all;
end