function [Fnulls, savedrand] = fdaComputeWnull(data,mList,permN,test)

if ~exist('permN','var')
    permN = 10000;
end

if ~exist('test','var')
    test = 'mixed';
end

if strcmp(test,'mixed')
    Fnulls = zeros(permN,size(data,1),3);       % we look for 2 main effects and 1 interaction
end

condN = size(data,3);
subjN = size(data,2);
savedrand = zeros(subjN,permN);
parfor perm = 1:permN
    fprintf('permutation %d/%d\n',perm,permN)
    locList = mList;
    if strcmp(test,'mixed')
        mListShuffled = locList(randperm(length(locList)));
        condShuffling = randi(condN, subjN, 1) == 2;
        savedrand(:,perm) = condShuffling;
        thisdata = data;
        thisdata(:,condShuffling,:) = thisdata(:,condShuffling,[2 1]);  % if more conditions make sure a random permutation
        [x, y, z] = fdaComputeFtest_mixedanova(thisdata,mListShuffled);
        Fnulls(perm,:,:) = [x; y; z]';
    end
end
