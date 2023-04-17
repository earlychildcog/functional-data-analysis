function [F_b, F_a] = fdaComputeF(data,condList)
% [F_b F_a] = computeF(data,condList)
% -----------------------------------
% compute functional F statistic for a linear model 
% Y = a + b*X + e
% where X is categorical (anova)
% 
% F_b:      F-statistic for slope
% F_a:      F-statistic for intercept
% data:     array of time x subject
% condList: array of subj x 1

F_b = zeros(1,size(data,1));
condNames = unique(condList);
condInWords = iscell(condList) && ischar(condList{1});
for t = 1:size(data,1)
    
    if condInWords

        thisdata{1} = data(t,strcmp(condList,condNames{1}));
        thisdata{2} = data(t,strcmp(condList,condNames{2}));
    else
        thisdata{1} = data(t,condList == condNames(1));
        thisdata{2} = data(t,condList == condNames(2));

    end
    
    K = length(condNames);
    
    
    grandMean = mean(data(t,:));                        % compute general mean
    groupMean = cellfun(@mean,thisdata);                % compute group means
    groupSize = cellfun(@length,thisdata);              % compute N for each group

    % compute sum of squares/degrees of freedom for the model
    SS_M = sum(groupSize.*(groupMean - grandMean).^2);  
    df_M = K - 1;

    % compute sum of squares/degrees of freedom for the error term
    SS_R = 0;
    df_R = 0;
    for k = 1:K
        df_R = df_R + length(thisdata{k}) - 1;
        for i = 1:length(thisdata{k})
            SS_R = SS_R + (thisdata{k}(i) - groupMean(k))^2;
        end
    end
    
    % check if we did things correct
    SS_T = sum((data(t,:) - grandMean).^2);
    df_T = length(data(t,:)) - 1;
    assert(abs(SS_T - (SS_R + SS_M)) < 10^(-6),'error in computation')
    assert(df_T == df_R + df_M,'error in computation')
    
    % compute mean squares and F ratio for main term
    MS_M = SS_M/df_M;
    MS_R = SS_R/df_R;
    F_b(t) = MS_M/MS_R;
    
    % compute mean squares and F ratio for intercept term
    




    
end