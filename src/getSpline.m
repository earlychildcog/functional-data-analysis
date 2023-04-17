function ys = getSpline(x, y, knots)
% gives cubic spline approximation
arguments
    x double
    y double
    knots double = augknt(x(1):((x(end)-x(1))/40):x(end),4);
end

if length(x) ~= length(y)
    error("sizes provided differ, x %d y %d", length(x), length(y))
end

sp = spap2(knots,4,x,y);

ys = fnval(sp,x);