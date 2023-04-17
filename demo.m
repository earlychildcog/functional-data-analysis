% script for demonstrating a functional data analysis example

% coming...
%% add paths
addpath(genpath("src"))

groupFile = 'data/groups.csv';
csvfolder = 'data/csv/';

%% construct splines
window = 0;
if window == 0
    duration_original = [6.5,17.5];
elseif window == 3
    duration_original = [14.5,17.5];
elseif window == -1
    duration_original = [5,6.5];
end
if mod(diff(duration_original),1) == 0
    dur_dilation = 5;
elseif mod(diff(duration_original),1) == 0.5
    dur_dilation = 6;
end
duration = duration_original*dur_dilation;
fs = 0.002;

splineN = duration(2)-duration(1)-1+4;
pupbasis = create_bspline_basis(duration, splineN, 4);
display(pupbasis);

figure
plot(pupbasis)


tvec = [duration_original(1):fs:duration_original(2)];
basismatrix = eval_basis(tvec*dur_dilation, pupbasis);
%% plot some derivative
figure
hold on
for s = 1:splineN
    plot(tvec,basismatrix(:,s))
end

%% load pupil data

look_into_subjs = [];       % plot which subjects

filenames = arrayfun(@(x)x.name, dir(fullfile(csvfolder, '*.csv')), 'UniformOutput', false);

M = readtable(groupFile);      % mirror condition file

mList = zeros(size(filenames));

subjN = length(filenames);
baseline = [2 5];               % set baseline
condNames = {'high','low'};
condN = length(condNames);
GoodTrials = zeros(subjN,2);
s = 0;

% loop through subjects to grab data
for f = 1:subjN
    filename = filenames{f};
    P = readtable(fullfile(csvfolder,filename));     
    P(P.trial < 3, :) = [];                  % remove fam trials if there
    P(strcmp(P.condition,'true'),:) = [];    % delete true belief trials
    P(mod(P.time,2) == 1, :) = [];           % remove orphan timepoints (prob not needed anymore)
    P.time = P.time/1000;                    % time to seconds
    
    if min(arrayfun(@(x)length(unique(P.trial(strcmp(P.condition,x)))),condNames)) == 0
        mList(f) = [];
    else
        s =  s + 1;
        % compute baseline for each trial
        Pbase = varfun(@mean,P(P.time >= baseline(1) & P.time <= baseline(2),:),'InputVariables','pupil','GroupingVariables','trial');
        Pbase.Properties.VariableNames{end} = 'basepupil';
        Pbase.GroupCount = [];
        P = join(P,Pbase,'Keys','trial');
        P.pupil = P.pupil -P.basepupil;         % baseline data

        % restrict time in IP
        P(P.time > max(tvec),:) = [];
        P(P.time < min(tvec),:) = [];
        trialList = unique(P.trial);            
        condList = arrayfun(@(x)P.condition(P.trial == x & P.time == tvec(1)),trialList);
        mList(s) = M.mirror(strcmp(M.session, P.session{1}));

        % if we wanted to plot some or all subjects
        if contains(filename,arrayfun(@num2str,look_into_subjs,'UniformOutput',false))
            figure
            title(filename)
            for c = 1:2
                subplot(2,1,c)
                hold on
                trialList = unique(P.trial(strcmp(P.condition, condNames{c})))';
                for tr = trialList
                    plot(P.time(P.trial == tr), P.pupil(P.trial == tr),'.')
                end
                legend(arrayfun(@num2str,trialList,'UniformOutput',false),'Location','eastoutside');
                title([condNames{c} ' ' filename]);
                ylim([-1 1])
                xlabel('time (s)')
                ylabel('pupil diameter (mm)')
            end
        end
        % average data per condition
        Pave = varfun(@nanmean,P,'InputVariables','pupil','GroupingVariables',{'condition','time'});
        Pave.Properties.VariableNames{end} = 'pupil';
        for c = 1:length(condNames)
            GoodTrials(s,c) = max(Pave.GroupCount(strcmp(Pave.condition,condNames{c})));
            if size(unique(Pave.GroupCount(strcmp(Pave.condition,'high'))),1) + size(unique(Pave.GroupCount(strcmp(Pave.condition,'low'))),1) > 2
                Pave(Pave.GroupCount == 1,:) = [];
            end
            pupmat{c}(:,s) = Pave.pupil(strcmp(Pave.condition,condNames{c}));
        end
    end
end

sesions = unique(P.session);

%% get spline approximation



% compute coefficients
simplemethod = true;            % we use the simple method
if simplemethod
    for c = 1:length(condNames)
        pupcoef{c} = basismatrix\pupmat{c};
    end
else
    [fdobj, df, gcv, coef, SSE, penmat, y2cMap] = smooth_basis(tvec*dur_dilation, pupmat, pupbasis);
    pupcoef = fdobj.coef;
end

% we use the coefficients to compute the approximation
subjN = s;
for c = 1:condN
    for s = 1:subjN
        SPave{c}(:,s) = zeros(size(tvec))';
        for splineNo = 1:splineN
            SPave{c}(:,s) = SPave{c}(:,s) + pupcoef{c}(splineNo,s)*basismatrix(:,splineNo);
        end
    end
end
%% other stuff
% plot
if true
disseq = zeros(1,subjN);
for s = 1:subjN
    if contains(filenames{s},arrayfun(@num2str,look_into_subjs,'UniformOutput',false))
        figure
        subplot(2,1,1)
        title(sprintf('subj %s',filenames{s}),'Interpreter','none')
        hold on
        for c =  1:condN
    %         plot(tvec, SPave{c}(:,s))
            plot(tvec, pupmat{c}(:,s))

        end
        ylim([-1 1])
        legend(condNames);
        subplot(2,1,2)
        plot(tvec, abs(pupmat{2}(:,s) - pupmat{1}(:,s)))
        sqdis = sqrt(sum((pupmat{2}(:,s) - pupmat{1}(:,s)).^2));
        disseq(s) = sqdis;
        title(sprintf('goodtrials (%d,%d) sqdis %f',GoodTrials(s,:), sqdis),'Interpreter','none')
        ylim([0 1]);
    end
%     legend(condNames);
end

figure, scatter(mean(GoodTrials,2)',disseq,'.')
hold on
% figure, scatter(min(GoodTrials,[],2)',disseq,'.')

[r p] = corr(mean(GoodTrials,2),disseq','rows','pairwise')
% [r p] = corr(min(GoodTrials,[],2),disseq','rows','pairwise')
mdl = fitlm(mean(GoodTrials,2),disseq','VarNames',{'mean amount good trials','curve deviation'})
figure
plot(mdl)
mdl = fitlm(min(GoodTrials,[],2),disseq','VarNames',{'min amount good trials','curve deviation'})
figure
plot(mdl)

badsubj = disseq > 15;
end
%% if we want to change to original data
if false
    SPave = pupmat;
end
%% do a simple T-test for each condition separately (between groups)

for c = 1:1
    pupave = SPave{c};
    tvec0 = tvec;
    ip = tvec > 14.5 & tvec < 17;
    tvec = tvec(ip);
    pupave = pupave(ip,:);
    if false
        pupave = pupmat{c};
        E = readtable('data\graph_data.csv');
        E(E.time < duration_original(1)*1000-2, :) = [];
        E = E(strcmp(E.condition, condNames{c}),:);
        pupave = reshape(E.pupil,[],50);
    end
    figure
    hold on
    plot(tvec, nanmean(pupave(:,logical(mList)),2),'LineWidth',2)
    plot(tvec, nanmean(pupave(:,~logical(mList)),2),'LineWidth',2)
    legend({'recognisers','non-recognisers'})
    xlabel('time (s)')
    ylabel('baseline corrected pupil diameter (mm)')
    title('PETRA (high demand condition)')
    grid on

    % compute null distribution

    permN = 10000;

    T = computeT(pupave,mList);
    Tnulls = computeTnull(pupave,mList,permN);
    plotFDtimecourse(T, Tnulls, tvec, 0.95, 'F', false);
    title(sprintf('Functional t-test, condition %s',condNames{c}));
    ax = gca;
    grid off

    [C_max cluster_range] = compute_cluster(T, 2, 2)

    C_nulls = compute_cluster(Tnulls, 2, 2)
    C_nulls_sort = sort(C_nulls);
    Cthress = C_nulls_sort(ceil(length(C_nulls_sort)*0.95))
Tmax = C_max
Fnulls = C_nulls;
Fnsort = sort(Fnulls);
mmax = length(Fnsort);
figure;
histogram(Fnulls,[0:100:(max(Fnulls)+500)],'Normalization','probability','HandleVisibility','off');
hold on
xline(Fnsort(floor(0.95*mmax)), 'r--', 'LineWidth',2);
xline(Tmax, 'color',[0 0.5 0.1], 'LineWidth',2);
legend('95% quantile','cluster sum-T-statistic')
xlabel('sum-T-statistic values')
title('empirical permutation distribution of %s T-statistics',perm_measure)
    tvec = tvec0;
    %

end
%% pooled t test (pooling conditions, between groups)

pupave = (SPave{1}+SPave{2})/2;
pupave(:,45) = nanmean((SPave{1}+SPave{2})/2,2);
permN = 100000;

T = computeT(pupave,mList);
Tnulls = computeTnull(pupave,mList,permN);
d = computeDfun(tvec,pupave,mList);
[f, Qnull, Qmax] = plotFDtimecourse(T, Tnulls, tvec, 0.95, 'T', false);
title(sprintf('Functional t-test, pooled fb conditions'));

saveresults = false;
if saveresults
    if window == -1
        thisname = 'fda_ttest_earlyeffect';
    else
        thisname = 'combined';
    end
    saveas(f,sprintf('results/fda_ttest_%s.fig',thisname))
    try
        exportgraphics(f,sprintf('results/fda_ttest_%s.jpg',thisname),'Resolution',300)
    catch
        warning('plot not saved, old version of matlab, old export version used')
    end
end

% find significant times

sigtimes = (Qnull' <= abs(T));
sigtimes = tvec(diff(sigtimes) ~= 0)
%% check variances for anova test
check_vareq = true;
if check_vareq      

    pupave = zeros([size(SPave{1}) condN]);
    for c = 1:condN
        pupave(:,:,c) = SPave{c};
    %     pupave(:,:,c) = pupmat{c};
    end
    bList = [mList]';
    
    badsbj = any(isnan(pupave),[1 3]);
    
    % delete bad subj
    pupave(:,badsbj,:) = [];
    bList(badsbj) = [];
    bNames = unique(bList);
    cNames = [1 2];

    mdata = arrayfun(@(x,y)pupave(:,bList == x,y),[bNames' bNames'] ,[cNames; cNames],'UniformOutput',false);
    
    vdata = cellfun(@(x)var(x, [], 2),mdata,'UniformOutput',false);

    figure
    hold on
    arrayfun(@(m,c,n,d)plot(tvec, vdata{m,c}./vdata{n,d}), [1 1 1 2], [1 1 2 1], [1 2 2 2], [2 1 2 2])
    plot(xlim, [4 4],'k--')
    plot(xlim, 1./[4 4],'k--')
    plot(xlim, [4 4]/2,'--','Color',[1 1 1]/2)
    plot(xlim, 2./[4 4],'--','Color',[1 1 1]/2)
    plot(xlim, [0 0],'k-')

    vdata_all = cellfun(@(x)var(x, [], [2 3]),mdata,'UniformOutput',false);
    v = vdata_all{1}./vdata_all{2};
    figure, plot(tvec, v)
    hold on
    plot(xlim, [4 4],'k--')
    plot(xlim, 1./[4 4],'k--')
    plot(xlim, [4 4]/2,'--','Color',[1 1 1]/2)
    plot(xlim, 2./[4 4],'--','Color',[1 1 1]/2)
    plot(xlim, [0 0],'k-')
    
   
    p = zeros(size(tvec));
    for c = 1:2
        for ind = 1:length(tvec)
    %     tpoint = 17;
        %     ind = find(tvec == tpoint);
            thispup = pupave(ind, :, c);
            thispup = thispup(:);
            p(ind) = vartestn(thispup, bList,'Display','off','TestType','LeveneAbsolute');
        end
        figure
        plot(tvec,p)
        hold on
        plot(xlim, 0.05*[1 1],'k--')
    end


    for ind = 1:length(tvec)
%     tpoint = 17;
    %     ind = find(tvec == tpoint);
        thispup = mean(pupave(ind, :, :),3);
        thispup = thispup(:);
        [p(ind) stats(ind)] = vartestn(thispup, bList,'Display','off','TestType','LeveneQuadratic');
    end
    figure
    plot(tvec,p)
    hold on
    plot(xlim, 0.05*[1 1],'k--')

    FL = [stats.fstat];
    figure
    plot(tvec,FL)
end


%% full anova

% prepare and cleanup data
condN = 2;
pupave = zeros([size(SPave{1}) condN]);
for c = 1:condN
    pupave(:,:,c) = SPave{c};
%     pupave(:,:,c) = pupmat{c};
end
bList = [mList]';

badsbj = any(isnan(pupave),[1 3]);          % find subjects with missing values

pupave(:,badsbj,:) = [];                    % delete bad subj
bList(badsbj) = [];                         % delete bad subj in grouping array



% get F statistic
[F_b, F_w, F_bw] = computeF_mixedanova(pupave,bList);

% plot the differences
figure
hold on
plot(tvec, mean(pupave(:,logical(mList)),2))
plot(tvec, mean(pupave(:,~logical(mList)),2))
legend({'recognisers','non-recognisers'})

figure
hold on
plot(tvec, mean(pupave(:,:,1),2))
plot(tvec, mean(pupave(:,:,2),2))
legend({'high','low'})

% compute null distribution
permN = 10000;
[Fnulls savedrand] = computeFnull(pupave,bList,permN);
% title(sprintf('Functional ANOVA: F values for Mirror (main effect)'))

%% plot the values
alpha = 0.95;
statname = 'F';
names = {'mirror','demand','interaction'};
clear f
clear Qnull
saveresults = false;

if saveresults
    save('results/fda_anova.mat','F_b', 'F_w', 'F_bw', 'Fnulls', 'savedrand', 'permN')
end

k = 1;
[f{k}, Qnull{k}, Qmax{k}] = plotFDtimecourse(F_b, Fnulls(:,:,k), tvec, alpha, statname, false);
title('Functional ANOVA: mirror main effect')

if saveresults
    saveas(f{k},sprintf('results/fda_anova_%s.fig',names{k}))
    try
        exportgraphics(f{k},sprintf('results/fda_anova_%s.jpg',names{k}),'Resolution',300)
    catch
        warning('plot not saved, old version of matlab, old export version used')
    end
end

k = 2;
[f{k}, Qnull{k}, Qmax{k}] = plotFDtimecourse(F_w, Fnulls(:,:,k), tvec, alpha, statname, false);
title('Functional ANOVA: demand main effect')

if saveresults
    saveas(f{k},sprintf('results/fda_anova_%s.fig',names{k}))
    try
        exportgraphics(f{k},sprintf('results/fda_anova_%s.jpg',names{k}),'Resolution',300)
    catch
        warning('plot not saved, old version of matlab, old export version used')
    end
end
k = 3;
[f{k}, Qnull{k}, Qmax{k}] = plotFDtimecourse(F_bw, Fnulls(:,:,k), tvec, alpha, statname, false);
title('Functional ANOVA: interaction mirror x demand')

if saveresults
    saveas(f{k},sprintf('results/fda_anova_%s.fig',names{k}))
    try
        exportgraphics(f{k},sprintf('results/fda_anova_%s.jpg',names{k}),'Resolution',300)
    catch
        warning('plot not saved, old version of matlab, old export version used')
    end
end

% find significant times

sigtimes = (Qnull{1}' <= F_b);
sigtimes = tvec(diff(sigtimes) ~= 0)