function [x, y, th] = get_parabolic_nodes(a,Q)
% calculates the position vectors for Q+1 nodes (Q segments) discretising a parabola
% y=ax^2.

x       = linspace(-1,1,10000);
y       = a.*x.^2;

ds      = sqrt((x(2:end) - x(1:end-1)).^2 + (y(2:end) - y(1:end-1)).^2);
s       = [0,cumsum(ds)];
delta   = max(s)/Q;
ss      = 0:delta:max(s);

yy      = spline(s,y,ss);
xx      = spline(s,x,ss);
x = xx/max(s); y = yy/max(s);

dy      = diff(yy);
dx      = diff(xx);
th      = atan(dy./dx); th = th';

%% testing

X = [x;y];
dX = [diff(x);diff(y)];
arc = 0;
for k = 1:length(dX)
    nX(k)   = norm(dX(:,k));
    arc     = arc+nX(k);
end
% arc

end % function

