function Tnulls = fdaComputeTtestNull(data,opts)
arguments
    data double
    opts.test string {mustBeMember(opts.test, ["one-sample" "two-sample"])} = "one-sample"
    opts.permN double = 1000
    opts.groupList double = []
    opts.nullMean double = 0
    opts.verbose logical = false
end

permN = opts.permN;
test = opts.test;
groupList = opts.groupList;
nullMean = opts.nullMean;
verbose = opts.verbose;

if test == "two-sample"
    assert(~isempty(groupList), "please provide a list for group membership in two-sample test")
end
    
if test == "one-sample" && size(data,2) == 2
    data = diff(data,1,2);
end
    
Tnulls = zeros(permN,size(data,1));
if test == "one-sample"
    parfor perm = 1:permN
        if verbose
            fprintf('permutation %d/%d\n',perm,permN)
        end
        Tnulls(perm,:) = fdaComputeTtest((randi(2,1,size(data,2))*2-3).*data,nullMean = nullMean);
    end
else
    parfor perm = 1:permN
        if verbose
            fprintf('permutation %d/%d\n',perm,permN)
        end
        thisList = groupList;
        thisList = thisList(randperm(length(thisList)));
        Tnulls(perm,:) = fdaComputeTtest(data,groupList = thisList, test = "two-sample");
    end
end
