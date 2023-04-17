function [F_b, F_w, F_bw] = fdaComputeFtest_mixedanova(data,bList)
% [F_a, F_b, F_ab] = computeF_mixedanova(data,condList)
% -----------------------------------
% compute functional F statistic for mixed effects two-way anova 
% data should have no missing values
% but otherwise design may be unbalanced between groups
% 
% F_b:      F-statistic for between factor
% F_w:      F-statistic for within factor
% F_ab:     F-statistic for interaction
% data:     array of time x subject x condition
% bList:    between subjects condition, array of subj x 1 (numeric)

%% prepare data
check_output     = false;

F_b = zeros(1,size(data,1));
F_w = zeros(1,size(data,1));
F_bw = zeros(1,size(data,1));

bNames = unique(bList);

B = length(bNames);

W = size(data,3);               % number of (w-s) conditions
N = size(data,2)*W;             % subjects x w-s conditions
s = size(data,2);
n = arrayfun(@(k)sum(bList == bNames(k)), 1:B);

%% compute vectorised

%% s1
data_split = arrayfun(@(k)data(:,bList == bNames(k),:), 1:B, 'UniformOutput',false);        % split the data into the b-s groups
grandMean = mean(data,[2 3]);                                                               % compute general sum at each timepoint

%% bettween groups fixed effect estimation
% this part computes the first part of the anova, ie the between subjects part
% it corresponds to the model EY ~ b0 + b*B + e

groupMeanB = cell2mat(cellfun(@(x)mean(x,[2 3]),data_split,'UniformOutput',false));         % compute mean for each b-s group
groupMeanSW = mean(data,3);                                                                 % compute mean for each subject across w-s conditions

% compute SoS/DoF for between groups fixed effect term
SS_b = sum(W*n.*(groupMeanB-grandMean).^2, 2);                                              
df_b = B - 1;                                                                               

% compute SoS/DoF for between groups error term
SS_bs = sum(W*(groupMeanSW - grandMean).^2, 2);                                                
SS_Rb = SS_bs - SS_b;
df_Rb = s - B;

% compute sum squares and F ratio for between term
MS_b = SS_b/df_b;
MS_Rb = SS_Rb/df_Rb;
F_b = MS_b./MS_Rb;

%% within groups estimation
% this part computes the second part of the anova, ie the within subjects part
% it corresponds to the model Y2-Y1 ~ b0 + b1*W + b2*Iwb + e

groupMeanW = mean(data,2);                % compute group Sums
groupMeanBW = cell2mat(cellfun(@(x)mean(x,2),data_split,'UniformOutput',false));                % compute group x condition means: T x B x W array

% compute SoS/DoF for between conditions factor
SS_w = sum(s*(groupMeanW - grandMean).^2, 3);
df_w = W - 1;

% compute SoS/DoF for the interaction term
%     SS_bw = sum(n'.*(groupMeanBW - grandMean).^2,'all');
SS_bw = sum(n.*(groupMeanBW).^2,[2 3]) - N*grandMean.^2 - SS_w - SS_b;
df_bw = (B-1)*(W-1);

% compute SoS/DoF for the error term between conditions
SS_T = sum((data - grandMean).^2,[2 3]);
SS_Rbw = SS_T - SS_bs - (SS_bw + SS_w);
df_Rbw = N-s-(B-1)*W;



% compute F ratio for within term
MS_w = SS_w/df_w;
MS_Rbw = SS_Rbw/df_Rbw;
F_w = MS_w./MS_Rbw;

% compute F ratio for interaction term
MS_bw = SS_bw/df_bw;
F_bw = MS_bw./MS_Rbw;

%% check if results compatible with matlab's repeated measures anova
if check_output
    tpoint = 4;
    thisdata = zeros(size(data,[2 3]));
    thisdata(:,:) = data(tpoint,:,:);        % get data for this timepoint
    [tbl,rm] = simple_mixed_anova(thisdata, bList');
    tbl.SumSq(2:end) - [SS_b(tpoint); SS_Rb(tpoint); SS_w(tpoint); SS_bw(tpoint); SS_Rbw(tpoint)]
    assert(abs(tbl.SumSq(2:end) - [SS_b(tpoint); SS_Rb(tpoint); SS_w(tpoint); SS_bw(tpoint); SS_Rbw(tpoint)]) < 10^(-8))
end



