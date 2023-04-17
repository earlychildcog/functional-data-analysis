function [Fnulls, savedrand] = fdaComputeFtest_null(data,groupList,opts)
arguments
    data double
    groupList double
    opts.nPerm (1,1) double = 10000
    opts.verbose logical = false
end

verbose = opts.verbose;
nPerm = opts.nPerm;
Fnulls = zeros(nPerm,size(data,1),3);       % we look for 2 main effects and 1 interaction

nCond = size(data,3);
nSubj = size(data,2);
savedrand = zeros(nSubj,nPerm);
parfor perm = 1:nPerm
    if verbose
        fprintf('permutation %d/%d\n',perm,nPerm)
    end
    locList = groupList;
    mListShuffled = locList(randperm(length(locList)));
    condShuffling = randi(nCond, nSubj, 1) == 2;
    savedrand(:,perm) = condShuffling;
    thisdata = data;
    thisdata(:,condShuffling,:) = thisdata(:,condShuffling,[2 1]);  % if more conditions make sure a random permutation
    [x, y, z] = fdaComputeFtest_mixedanova(thisdata,mListShuffled);
    Fnulls(perm,:,:) = [x, y, z];
end
