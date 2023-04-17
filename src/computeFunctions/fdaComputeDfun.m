function d = fdaComputeDfun(tvec,data,condList)
d = zeros(size(tvec));
condNames = unique(condList);
condInWords = iscell(condList) && ischar(condList{1});
for t = 1:length(tvec)
    if condInWords

        data1 = data(t,strcmp(condList,condNames{1}));
        data2 = data(t,strcmp(condList,condNames{2}));
    else
        data1 = data(t,condList == condNames(1));
        data2 = data(t,condList == condNames(2));

    end
    
    dnom = mean(data1) - mean(data2);
    sv1      = sum((data1-mean(data1)).^2);
    sv2      = sum((data2-mean(data2)).^2);
    pooledSD =  sqrt((sv1 + sv2)/(length(data1) + length(data2)-2)); % pooled Standard Deviation
    ddenom   = pooledSD;

    
    d(t) = dnom/ddenom;
end