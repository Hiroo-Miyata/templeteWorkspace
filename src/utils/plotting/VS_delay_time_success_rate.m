function VS_delay_time_success_rate(failedLabels, rewardLabels, difficultyLabels, delayTimes, options)

    arguments
        failedLabels
        rewardLabels
        difficultyLabels
        delayTimes
        options.Label
        options.OutputFolder
        options.Ylim = [20 100]
    end

    nrewards = 3; rewards = [1 2 3];
    ndifficulties = 2; difficulties = [0 1];
    rewColors = [1 0 0; 1 0.6470 0; 0 0 1]; diffColors = [0 0.447 0.741; 0.466 0.674 0.188];
    direColors = {[1 .5 .5],[.75 .75 .5],[.5 1 .5],[.25 .75 .5],[0 .5 .5],[0.25 0.25 .5],[0.5 0 .5],[0.75 0.25 .5]};
    DiffStyle = ["-", ":"];
    DelayTimes = 300:50:1050; nDelayTimes = length(DelayTimes);

    figure; hold on;

    for j=1:ndifficulties
        for i=1:nrewards
            Y = zeros(nDelayTimes, 1);
            for k=1:nDelayTimes
                curInds = rewardLabels == rewards(i) & difficultyLabels == difficulties(j) & delayTimes == DelayTimes(k);
                curInds2 = curInds & failedLabels;
                Y(k) = 100 * (1-sum(curInds2)/sum(curInds));
            end
            h(j) = plot(Y, 'Color', rewColors(i, :), 'LineWidth', 2.5,  'LineStyle', DiffStyle(j));
        end
    end
    set(gca, 'fontsize', 16, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
    xlim([0.7 nDelayTimes+.3]); ylim(options.Ylim); xticks(1:3:nDelayTimes); xticklabels(DelayTimes(1:3:end));
    xlabel("Delay Length"); ylabel(options.Label); legend(h, ["Tiny", "Huge"], Location="best"); set(gcf,'position',[0,0,550,550]);
    saveas(gcf, options.OutputFolder+"-vs-DelayLength.jpg");
    close all;
end