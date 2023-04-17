function [f, Qnull, QnullMaxLim] = fdaPlotTimecourse(F, Fnulls, tvec, opts)
arguments
    F
    Fnulls
    tvec            double  = 1:length(F)
    opts.alpha      (1,1) double = 0.95
    opts.statistic   char {mustBeMember(opts.statistic, {'F' 'T'})} = 'T'
    opts.control    string {mustBeMember(opts.control, ["point" "max" "theoretical"])} = "point"
    opts.sided      double = 2
    opts.theorThres double = [];
end
% [f, Qnull, QnullMaxLim] = fdaPlotTimecourse(F, Fnulls, tvec, opts)
% generic function to plot any functional statistic
% dispmax = 0 (default) to not use max distribution threshold in plot
% sided = 1 for one sided (default), 2 for two-sided (alpha remains the same).


alpha = opts.alpha;
statistic = opts.statistic;
control = opts.control;
sided = opts.sided;
theorThres = opts.theorThres;
if isempty(theorThres)
        if statistic == 'F'
            theorThres = finv(alpha, 1,Inf);
        elseif statistic == 'T'
            theorThres = tinv(alpha, Inf);
        end
end




f = figure;
hold on
plot(tvec,F,'LineWidth',2)
if sided == 1
    Qnull = quantile(Fnulls,alpha);
    leg = [sprintf("%s-statistic function",statistic)];
    if any(control == "point")
        plot(tvec,Qnull,'LineWidth',2)
        leg = [leg sprintf("upper %.3f quantile null distribution",(1-alpha))];
    end
    FnullDiff = max(Fnulls,[],2);
    QnullMax = quantile(FnullDiff,alpha);
    QnullMaxLim = QnullMax;

    oldylim = ylim;
    ylim(oldylim + [0 0.2].*oldylim)
    if any(control == "max")
        % plot maximal thresholds
        plot(xlim,[QnullMax QnullMax],'k--')
        leg = [leg sprintf("max null distribution for alpha %.3f",alpha)];
    end
    if any(control == "theoretical")
        plot(xlim,[theorThres theorThres],'k--')
        leg = [leg sprintf("theoretical threshold for alpha %.3f",alpha)];
    end
        
    
    legend(leg,'Location','best')
elseif sided == 2
    % plot pointwise thresholds
    QnullPos = quantile(Fnulls,1-(1-alpha)/2);
    QnullNeg = quantile(Fnulls,(1-alpha)/2);
    leg = [sprintf("%s-statistic function",statistic)];
    if any(control == "point")
        plot(tvec,QnullPos,'r','LineWidth',2)
        plot(tvec,QnullNeg,'r','LineWidth',2)
        leg = [leg sprintf("%.3f/%.3f quantile null distribution",(1-alpha)/2,1-(1-alpha)/2)];
    end
    FnullPos = max(Fnulls,[],2);
    QnullMax = quantile(FnullPos,1-(1-alpha)/2);
    FnullNeg = min(Fnulls,[],2);
    QnullMin = quantile(FnullNeg,(1-alpha)/2);

    QnullMaxLim = [QnullMin QnullMax];

    oldylim = ylim;
    ylim(oldylim + [0 0.2].*oldylim)
    if any(control == "max")
        % plot maximal thresholds
        plot(xlim,[QnullMax QnullMax],'k--')
        plot(xlim,[QnullMin QnullMin],'k--')
        leg = [leg sprintf("limits of min/max null distribution for alpha %.3f",alpha)];
    end
    if any(control == "theoretical")
        plot(xlim,[theorThres theorThres],'k--')
        plot(xlim,-[theorThres theorThres],'k--')
        leg = [leg sprintf("theoretical threshold for alpha  %.3f",alpha)];
    end
        
    
    legend(leg,'Location','best')
    
else
    error('sided must be 1 or 2')
end

xlabel('time (s)')
ylabel(sprintf('%s-values',statistic))    
xrange = xlim;
xticks(xrange(1):1:xrange(end))
grid on
title(sprintf('Functional %s-test',statistic))



if any(control == "point")        % to add shading later
    if sided == 1
        X = find(F >= Qnull');
    elseif sided == 2
        X = find(F >= QnullPos' | F<=QnullNeg');
    end
elseif any(control == "max")
    if sided == 1
        X = find(F >= QnullMax);
    elseif sided == 2
        X = find(F >= QnullMax | F <= QnullMin);
    end
elseif any(control == "theoretical")
    if sided == 1
        X = find(F >= theorThres);
    elseif sided == 2
        X = find(F >= theorThres | F<=-theorThres);
    end
end
yrange = ylim;
for s = 1:length(X)
    I = [X(s) X(s)+1];
    fill(([I flip(I)]*2)/1000 + tvec(1), [yrange(1) yrange(1) yrange(2) yrange(2)],'r', 'facealpha', 0.2, 'linestyle','none')
end

legend(leg, 'Location','northeast')

xlim([min(tvec) max(tvec)]);
grid off
end