function T = fdaComputeTtest(data,opts)
% compute functional T statistic
% data is number array time X 1
% option can be:
%       1. numeric array of two variables for two-sample t-test
%       2. a real number for the null hypothesis of one-sample t-test (or empty for mean zero)
arguments
    data double
    opts.groupList double = []
    opts.test string {mustBeMember(opts.test, ["one-sample" "two-sample"])} = "one-sample"
    opts.nullMean double = 0
    opts.groupNames double = []
end

test = opts.test;
groupList = opts.groupList;
nullMean = opts.nullMean;
groupNames = opts.groupNames;
if isempty(groupNames) && test == "two-sample"
    groupNames = unique(groupList);
end


T = zeros(size(data,1),1);

% in case of one sample t-test check if data is one sample, else take difference
if test == "one-sample" && size(data,3) == 2
    data = diff(data,1,3);
end
if test == "two-sample"
    
    
    % perform two-sample t-test
    data1 = data(:,groupList == groupNames(1));
    data2 = data(:,groupList == groupNames(2));
    n1    = sum(~isnan(data1),2);
    n2    = sum(~isnan(data2),2);
    Tnom = (mean(data1,2,"omitnan") - mean(data2,2,"omitnan"));
    Tdenom = sqrt(var(data1,[],2,"omitnan")./n1 + var(data2,[],2,"omitnan")./n2);
    T = Tnom./Tdenom;
else

    % perform a one-sample t-test
    
    Tnom = mean(data,2) - nullMean;
    Tdenom = std(data,[],2)/sqrt(size(data,2));
    T = Tnom./Tdenom;
        
end