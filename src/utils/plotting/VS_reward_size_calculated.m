function VS_reward_size_calculated(Y, Yerr, options)

    arguments
        Y (3, 2) double
        Yerr = nan
        options.Label (1, :) char
        options.OutputFolder (1, :) char
    end

    label = options.Label; outputFolder = options.OutputFolder;

    nrewards = 3; rewards = [1 2 3];
    ndifficulties = 2; difficulties = [0 1];
    rewColors = [1 0 0; 1 0.6470 0; 0 0 1]; diffColors = [0 0.447 0.741; 0.466 0.674 0.188];
    direColors = {[1 .5 .5],[.75 .75 .5],[.5 1 .5],[.25 .75 .5],[0 .5 .5],[0.25 0.25 .5],[0.5 0 .5],[0.75 0.25 .5]};
    DiffStyle = ["-", ":"];

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