function [Qnull ] = fdaPlotNullDist(S, Snull, opts)
arguments
    S (1,1) double
    Snull double
    opts.type string {mustBeMember(opts.type, ["max" "cluster" "average"])} = "cluster"
    opts.transformation string {mustBeMember(opts.transformation, ["none", "log"])} = "none"
    opts.alpha double = 0.95
end
if opts.transformation == "log"
    S = log10(S + 1);
    Snull = log10(Snull + 1);
    transftitl = "(log tranformed)";
    scale = " (log scale)";
else
    transftitl = "";
    scale = "";
end


nullThres = quantile(Snull, opts.alpha);
step = (max(Snull) - min(Snull))/100;


figure;
histogram(Snull,[0:step:max(Snull)],'Normalization','probability','HandleVisibility','off');
hold on
xline(nullThres, 'r--', 'LineWidth',2);
xline(S, 'color',[0 0.5 0.1], 'LineWidth',2);
legend(sprintf("%.2f quantile threshold",opts.alpha), opts.type + " statistic")
xlabel(opts.type + " statistic values" + scale)
title("empirical null " + opts.type + " distribution " + transftitl)


if opts.transformation == "log"
    xx = xlim;
    xticks(log10([1, 10.^[0:floor(xx(2))]+1]))
    xticklabels(num2cell(round(10.^get(gca,'XTick')-1)));
end