function [C_max cluster_range] = fdaComputeCluster(F, opts)
arguments
    F
    opts.thres (1,1) double = 1.9600
    opts.sided (1,1) {mustBeMember(opts.sided, [-1, 1, 2])} = 2
    opts.type string {mustBeMember(opts.type, ["t-sum" "cluster-extent"])} = "t-sum"
end
% function [C_sum cluster_range] = compute_cluster(F, thres, sided)
%
% computes the highest cluster above given threshold
%
% Output arguments:
% C_sum:            cluster statistic
% cluster_range:    if data is one dimensional, gives the extent of the detected cluster
%
% Input arguments:
% F:                the statistic data in row data; if 2-dimensional, computes a cluster statistic for each column 
% thress:           the forming clusters thresshold (should always be positive)
% sided:            2 for looking both positive and negative clusters. 1 for positive only and -1 for negative only.

%% input arguments processing
thres = opts.thres;
sided = opts.sided;
switch opts.type
    case "t-sum"
        normExponent = 1;
    case "cluster-extent"
        normExponent = 0;
end

% if two sided then we care about both negative and positive clusters; hence, due to continuity, we can take absolute values
if sided == 2
    F = abs(F);
elseif sided == -1
    F = -F;
end

% should change in future
if size(F,2) > 1
    F = F';
end

%% get clusters

nPerm = size(F,2);
F = [zeros(1,nPerm); F; zeros(1,nPerm)];                  % in case a cluster starts in the beginning or ends in the end of the time period
C = [F > thres];
C_diff = [diff(C); zeros(1,nPerm)];
C_start = C_diff == 1;
C_stat = (F.*C);
if normExponent > 1
    C_stat = C_stat.^normExponent;
elseif normExponent == 0
    C_stat = double(C_stat ~= 0);
end

C_sum = cumsum(C_stat);
C_sum_norm = C_sum - cumsum(C_sum.*C_start);
C_max = max(C_sum_norm);
        

if size(F,2) > 1
    cluster_range = [];
else
    Cstart_num = find(C_start);
    if ~isempty(Cstart_num)
        cluster_end = find(C_sum_norm.*C == C_max);
        cluster_start = Cstart_num(Cstart_num < cluster_end);
        cluster_start = cluster_start(end);
        cluster_range = [cluster_start cluster_end];
    else
        cluster_range = [NaN NaN];
    end
end




