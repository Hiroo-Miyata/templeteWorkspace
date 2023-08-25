function VS_direction_success_rate(failedLabels, rewardLabels, difficultyLabels, directionLabels, options)

    arguments
        failedLabels
        rewardLabels
        difficultyLabels
        directionLabels
        options.Label
        options.OutputFolder
        options.Ylim = [20 100]
        options.DiffLegendPos = "southwest"
        options.RewLegendPos = "south"
    end


    [rewardNames, rewardLegends, rewColors, diffColors, direColors, DiffStyle, DelayTimes, nDelayTimes] = getExperimentConstants();
    rewards = unique(rewardLabels); nrewards = length(rewards);
    difficulties = unique(difficultyLabels); ndifficulties = length(difficulties);
    directions = unique(directionLabels); ndirections = length(directions);

    figure; hold on;
    for j=1:ndifficulties
        for i=1:nrewards
            Y = zeros(ndirections, 1);
            for k=1:ndirections
                curInds = rewardLabels == rewards(i) & difficultyLabels == difficulties(j) & directionLabels == directions(k);
                curInds2 = curInds & failedLabels;
                Y(k) = 100 * (1-sum(curInds2)/sum(curInds));
            end
            lh(j) = plot(Y, 'Color', rewColors(rewards(i), :), 'LineWidth', 2.5,  'LineStyle', DiffStyle(j));
            lr(i) = lh(j);
        end
    end
    set(gca, 'fontsize', 16, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
    xlim([0.7 ndirections+.3]); ylim(options.Ylim); xticks(1:ndirections); xticklabels([0 45 90 135 180 225 270 315]);
    xlabel("Direction"); ylabel(options.Label); set(gcf,'position',[0,0,550,550]);
    legend(lh, ["Tiny", "Huge"], Location=options.DiffLegendPos); ah1=axes('position',get(gca,'position'),'visible','off');
    leg2=legend(ah1,lr,rewardLegends, Location=options.RewLegendPos);
    set(gca, 'fontsize', 20, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
    saveas(gcf, options.OutputFolder+"-vs-Direction.jpg");
    close all;
end