function corrHistogram(Y, options)
% exmaple
% corrHistogram(Yall(:), "Label", "Correlation Coefficient", "Title", "RT vs TriggerAxis-at-GC", "OutputFolder", outputFolder+"all-at-GC")
    arguments
        Y (:, 1) double
        options.Label (1,:) char
        options.OutputFolder (1,:) char
        options.Title (1,:) char = nan
        options.Xlim (1,:) char {mustBeMember(options.Xlim, ["none", "corrcoeff"])} = "corrcoeff"
    end

    figure; hold on;
    histogram(Y, 'FaceColor', 'k');
    
    p = signrank(Y(:));
    % add median line with red color
    line([median(Y) median(Y)], ylim, 'Color', 'r', 'LineWidth', 2);
    % show the median correlation coefficient, but is rounded to 3 decimal
    % show the p value with a * 10^b format, but a is rounded to 2 decimal
    b = floor(log10(p));
    a = round(p/(10^b), 2);
    if median(Y) < 0.75 * (max(xlim) - min(xlim)) + min(xlim)
        text(median(Y), 0.8*max(ylim), " Median: " + num2str(round(median(Y), 3)) + " ", 'Color', 'r', 'fontsize', 24);
        text(median(Y), 0.7*max(ylim), " p: " + num2str(a) + " * 10^" + num2str(b) + " ", 'Color', 'r', 'fontsize', 24);
    else
        % write left side
        text(median(Y), 0.8*max(ylim), " Median: " + num2str(round(median(Y), 3)) + " ", 'Color', 'r', 'fontsize', 24, 'HorizontalAlignment', 'right');
        text(median(Y), 0.7*max(ylim), "p: " + num2str(a) + " * 10^" + num2str(b), 'Color', 'r', 'fontsize', 24, 'HorizontalAlignment', 'right');
    end
    % show the black dot line of correlation coefficient = 0
    line([0 0], ylim, 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--');
    xlabel(options.Label);
    ylabel("Counts");
    if options.Xlim == "corrcoeff"
        xlim([-1 1]);
    end
    if ~isnan(options.Title)
        title(options.Title);
    end
    set(gca, 'fontsize', 20, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
    set(gcf,'position',[0,0,550,550]);
    saveas(gcf, options.OutputFolder+"-correlation-histogram.png");
    close all
end