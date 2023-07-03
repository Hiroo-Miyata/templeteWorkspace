function VS_reward_size_success_rate(failedLabels, rewardLabels, difficultyLabels, options)

    arguments
        failedLabels
        rewardLabels
        difficultyLabels
        options.Label
        options.OutputFolder
        options.Ylim = [20 100]
    end


    nrewards = 3; rewards = [1 2 3];
    ndifficulties = 2; difficulties = [0 1];
    rewColors = [1 0 0; 1 0.6470 0; 0 0 1]; diffColors = [0 0.447 0.741; 0.466 0.674 0.188];
    direColors = {[1 .5 .5],[.75 .75 .5],[.5 1 .5],[.25 .75 .5],[0 .5 .5],[0.25 0.25 .5],[0.5 0 .5],[0.75 0.25 .5]};
    DiffStyle = ["-", ":"];

    figure; hold on;
    for j=1:ndifficulties
        Y = zeros(3,1);
        for i=1:nrewards
            curInds = rewardLabels == rewards(i) & difficultyLabels == difficulties(j);
            curInds2 = curInds & failedLabels;
            Y(i) = 100 * (1-sum(curInds2)/sum(curInds));
        end
        plot(Y, 'Color', diffColors(j, :), 'LineWidth', 2.5);
    end
    set(gca, 'fontsize', 16, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
    xlim([0.7 3.3]); ylim(options.Ylim); xticks(1:3); xticklabels(["S", "M", "L"]);
    ylabel(options.Label); legend(["Tiny", "Huge"], Location="best"); set(gcf,'position',[0,0,550,550]);
    saveas(gcf, options.OutputFolder+"-vs-Reward.jpg"); close all;
end