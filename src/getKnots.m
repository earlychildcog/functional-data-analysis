function knots = getKnots(x, splineN)
% get knots by number of splines


xstep = (x(end) - x(1))/(splineN-3);
xvec = x(1):xstep:x(end);
knots = augknt(xvec,4);