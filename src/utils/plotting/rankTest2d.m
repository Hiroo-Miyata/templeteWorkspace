function rankTest2d(x, y, options)

    arguments
        x (:,1) double
        y (:,1) double
        options.XLabel (1,:) char
        options.YLabel (1,:) char
        options.AxisColor (1,:) char {mustBeMember(options.AxisColor, ["none", "diffColor"])} = "diffColor"
        options.OutputFolder (1,:) char
    end

    rewColors = [1 0 0; 1 0.6470 0; 0 0 1]; diffColors = [0 0.447 0.741; 0.466 0.674 0.188];

    figure; hold on;
    for i = 1:length(x)
        scatter(x, y, 150, 'k', "marker", "*");
    end
    xy = cat(1, x, y);
    maxlim = max(xy) + 0.1 *(max(xy) - min(xy));
    minlim = min(xy) - 0.1 *(max(xy) - min(xy));
    plot([minlim maxlim], [minlim maxlim], 'k', 'LineWidth', 2);
    meanX = mean(x); meanY = mean(y);
    scatter(meanX, meanY, 450, 'r', 'marker', 'pentagram', 'LineWidth', 2);
    % project the mean point onto the line x = y by red line
    dot_product = dot([1 1], [meanX meanY]) / 2;
    plot([dot_product meanX], [dot_product meanY], 'r', 'LineWidth', 3);

    diffxy = y - x;
    p = signrank(diffxy(:));
    if p < 0.001
        textp = "***";
    elseif p < 0.01
        textp = "**";
    elseif p < 0.05
        textp = "*";
    else
        textp = "n.s.";
    end
    text(mean([dot_product meanX]), mean([dot_product meanY])*1.07, textp, 'Color', 'r', 'FontSize', 40);

    lx = xlabel(options.XLabel);
    ly = ylabel(options.YLabel);
    if options.AxisColor == "diffColor"
        lx.Color = diffColors(1,:);
        ly.Color = diffColors(2,:);
    end
    set(gca, 'fontsize', 20, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
    set(gcf,'position',[0,0,550,550]); xlim([minlim maxlim]); ylim([minlim maxlim]);
    saveas(gcf, options.OutputFolder+"-signTest-2d.jpg");
    close all
end