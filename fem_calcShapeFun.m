% Calculates shape functions and derivatives for both velocity and
% pressure on [-1,1]x[-1,1]x[-1,1] and evaluates them at the Gaussian
% quadrature nodes given in GQ
%
% [shape, dshape] = fem_calcShapeFun(GQ)
%
% ========================================================================
% Inputs
% ========================================================================
% GQ = structure containing required Gaussian quadrature nodes and weights
%
% ========================================================================
% Outputs
% ========================================================================
% shape.v = velocity shape functions at GQ nodes
% shape.p = pressure shape functions at GQ nodes
% dshape.v = derivatives of velocity shape functions at GQ nodes
% dshape.p = derivatives of pressure shape functions at GQ nodes

function [shape, dshape] = fem_calcShapeFun(GQ)

for k=1:8
    s = GQ.node(k,1);
    t = GQ.node(k,2);
    r = GQ.node(k,3);
    
    % velocity shape functions
    shape.v(1,k) = 1/8*s*(s-1)*t*(t-1)*r*(r-1);
    shape.v(7,k) = 1/8*s*(s+1)*t*(t-1)*r*(r-1);
    shape.v(9,k) = 1/8*s*(s+1)*t*(t+1)*r*(r-1);
    shape.v(3,k) = 1/8*s*(s-1)*t*(t+1)*r*(r-1);
    shape.v(19,k) = 1/8*s*(s-1)*t*(t-1)*r*(r+1);
    shape.v(25,k) = 1/8*s*(s+1)*t*(t-1)*r*(r+1);
    shape.v(27,k) = 1/8*s*(s+1)*t*(t+1)*r*(r+1);
    shape.v(21,k) = 1/8*s*(s-1)*t*(t+1)*r*(r+1);
    shape.v(4,k) = 1/4*(1-s^2)*t*(t-1)*r*(r-1);
    shape.v(8,k) = 1/4*s*(s+1)*(1-t^2)*r*(r-1);
    shape.v(6,k) = 1/4*(1-s^2)*t*(t+1)*r*(r-1);
    shape.v(2,k) = 1/4*s*(s-1)*(1-t^2)*r*(r-1);
    shape.v(10,k) = 1/4*s*(s-1)*t*(t-1)*(1-r^2);
    shape.v(16,k) = 1/4*s*(s+1)*t*(t-1)*(1-r^2);
    shape.v(18,k) = 1/4*s*(s+1)*t*(t+1)*(1-r^2);
    shape.v(12,k) = 1/4*s*(s-1)*t*(t+1)*(1-r^2);
    shape.v(22,k) = 1/4*(1-s^2)*t*(t-1)*r*(r+1);
    shape.v(20,k) = 1/4*s*(s-1)*(1-t^2)*r*(r+1);
    shape.v(24,k) = 1/4*(1-s^2)*t*(t+1)*r*(r+1);
    shape.v(26,k) = 1/4*s*(s+1)*(1-t^2)*r*(r+1);
    shape.v(5,k) = 1/2*(1-s^2)*(1-t^2)*r*(r-1);
    shape.v(13,k) = 1/2*(1-s^2)*t*(t-1)*(1-r^2);
    shape.v(17,k) = 1/2*s*(s+1)*(1-t^2)*(1-r^2);
    shape.v(15,k) = 1/2*(1-s^2)*t*(t+1)*(1-r^2);
    shape.v(11,k) = 1/2*s*(s-1)*(1-t^2)*(1-r^2);
    shape.v(23,k) = 1/2*(1-s^2)*(1-t^2)*r*(r+1);
    shape.v(14,k) = (1-s^2)*(1-t^2)*(1-r^2);
    
    % derivatives of velocity shape functions w.r.t s
    dshape.v(1,1,k) = 1/8*(2*s-1)*(t-1)*t*(r-1)*r;
    dshape.v(1,7,k) = 1/8*(2*s+1)*(t-1)*t*(r-1)*r;
    dshape.v(1,9,k) = 1/8*(2*s+1)*t*(t+1)*(r-1)*r;
    dshape.v(1,3,k) = 1/8*(2*s-1)*t*(t+1)*(r-1)*r;
    dshape.v(1,19,k) = 1/8*(2*s-1)*t*(t-1)*(r+1)*r;
    dshape.v(1,25,k) = 1/8*(2*s+1)*t*(t-1)*(r+1)*r;
    dshape.v(1,27,k) = 1/8*(2*s+1)*t*(t+1)*(r+1)*r;
    dshape.v(1,21,k) = 1/8*(2*s-1)*t*(t+1)*(r+1)*r;
    dshape.v(1,4,k) = (-1)*1/2*s*(t-1)*t*(r-1)*r;
    dshape.v(1,8,k) = (-1)*1/4*(2*s+1)*(t^2-1)*(r-1)*r;
    dshape.v(1,6,k) = (-1)*1/2*s*t*(t+1)*(r-1)*r;
    dshape.v(1,2,k) = (-1)*1/4*(2*s-1)*(t^2-1)*(r-1)*r;
    dshape.v(1,10,k) = (-1)*1/4*(2*s-1)*(t-1)*t*(r^2-1);
    dshape.v(1,16,k) = (-1)*1/4*(2*s+1)*(t-1)*t*(r^2-1);
    dshape.v(1,18,k) = (-1)*1/4*(2*s+1)*(t+1)*t*(r^2-1);
    dshape.v(1,12,k) = (-1)*1/4*(2*s-1)*(t+1)*t*(r^2-1);
    dshape.v(1,22,k) = (-1)*1/2*s*(t-1)*t*r*(r+1);
    dshape.v(1,20,k) = (-1)*1/4*(2*s-1)*(t^2-1)*r*(r+1);
    dshape.v(1,24,k) = (-1)*1/2*s*t*(t+1)*r*(r+1);
    dshape.v(1,26,k) = (-1)*1/4*(2*s+1)*(t^2-1)*r*(r+1);
    dshape.v(1,5,k) = s*(t^2-1)*(r-1)*r;
    dshape.v(1,13,k) = s*(t-1)*t*(r^2-1);
    dshape.v(1,17,k) = 1/2*(2*s+1)*(t^2-1)*(r^2-1);
    dshape.v(1,15,k) = s*t*(t+1)*(r^2-1);
    dshape.v(1,11,k) = 1/2*(2*s-1)*(t^2-1)*(r^2-1);
    dshape.v(1,23,k) = s*(t^2-1)*r*(r+1);
    dshape.v(1,14,k) = (-1)*2*s*(t^2-1)*(r^2-1);
    
    % derivatives of velocity shape functions w.r.t t
    dshape.v(2,1,k) = 1/8*(s-1)*s*(2*t-1)*(r-1)*r;
    dshape.v(2,7,k) = 1/8*s*(s+1)*(2*t-1)*(r-1)*r;
    dshape.v(2,9,k) = 1/8*s*(s+1)*(2*t+1)*(r-1)*r;
    dshape.v(2,3,k) = 1/8*s*(s-1)*(2*t+1)*(r-1)*r;
    dshape.v(2,19,k) = 1/8*s*(s-1)*(2*t-1)*(r+1)*r;
    dshape.v(2,25,k) = 1/8*s*(s+1)*(2*t-1)*(r+1)*r;
    dshape.v(2,27,k) = 1/8*s*(s+1)*(2*t+1)*(r+1)*r;
    dshape.v(2,21,k) = 1/8*s*(s-1)*(2*t+1)*(r+1)*r;
    dshape.v(2,4,k) = (-1)*1/4*(s^2-1)*(2*t-1)*(r-1)*r;
    dshape.v(2,8,k) = (-1)*1/2*s*(s+1)*t*(r-1)*r;
    dshape.v(2,6,k) = (-1)*1/4*(s^2-1)*(2*t+1)*(r-1)*r;
    dshape.v(2,2,k) = (-1)*1/2*(s-1)*s*t*(r-1)*r;
    dshape.v(2,10,k) = (-1)*1/4*(s-1)*s*(2*t-1)*(r^2-1);
    dshape.v(2,16,k) = (-1)*1/4*s*(s+1)*(2*t-1)*(r^2-1);
    dshape.v(2,18,k) = (-1)*1/4*s*(s+1)*(2*t+1)*(r^2-1);
    dshape.v(2,12,k) = (-1)*1/4*s*(s-1)*(2*t+1)*(r^2-1);
    dshape.v(2,22,k) = (-1)*1/4*(s^2-1)*(2*t-1)*r*(r+1);
    dshape.v(2,20,k) = (-1)*1/2*s*(s-1)*t*r*(r+1);
    dshape.v(2,24,k) = (-1)*1/4*(s^2-1)*(2*t+1)*r*(r+1);
    dshape.v(2,26,k) = (-1)*1/2*s*(s+1)*t*r*(r+1);
    dshape.v(2,5,k) = (s^2-1)*t*(r-1)*r;
    dshape.v(2,13,k) = 1/2*(s^2-1)*(2*t-1)*(r^2-1);
    dshape.v(2,17,k) = s*(s+1)*t*(r^2-1);
    dshape.v(2,15,k) = 1/2*(s^2-1)*(2*t+1)*(r^2-1);
    dshape.v(2,11,k) = (s-1)*s*t*(r^2-1);
    dshape.v(2,23,k) = (s^2-1)*t*r*(r+1);
    dshape.v(2,14,k) = (-2)*(s^2-1)*t*(r^2-1);
    
    % derivatives of velocity shape functions w.r.t r
    dshape.v(3,1,k) = 1/8*(s-1)*s*(t-1)*t*(2*r-1);
    dshape.v(3,7,k) = 1/8*s*(s+1)*(t-1)*t*(2*r-1);
    dshape.v(3,9,k) = 1/8*s*(s+1)*t*(t+1)*(2*r-1);
    dshape.v(3,3,k) = 1/8*s*(s-1)*t*(t+1)*(2*r-1);
    dshape.v(3,19,k) = 1/8*s*(s-1)*t*(t-1)*(2*r+1);
    dshape.v(3,25,k) = 1/8*s*(s+1)*t*(t-1)*(2*r+1);
    dshape.v(3,27,k) = 1/8*s*(s+1)*t*(t+1)*(2*r+1);
    dshape.v(3,21,k) = 1/8*s*(s-1)*t*(t+1)*(2*r+1);
    dshape.v(3,4,k) = (-1)*1/4*(s^2-1)*(t-1)*t*(2*r-1);
    dshape.v(3,8,k) = (-1)*1/4*s*(s+1)*(t^2-1)*(2*r-1);
    dshape.v(3,6,k) = (-1)*1/4*(s^2-1)*t*(t+1)*(2*r-1);
    dshape.v(3,2,k) = (-1)*1/4*(s-1)*s*(t^2-1)*(2*r-1);
    dshape.v(3,10,k) = (-1)*1/2*(s-1)*s*(t-1)*t*r;
    dshape.v(3,16,k) = (-1)*1/2*s*(s+1)*(t-1)*t*r;
    dshape.v(3,18,k) = (-1)*1/2*s*(s+1)*(t+1)*t*r;
    dshape.v(3,12,k) = (-1)*1/2*s*(s-1)*(t+1)*t*r;
    dshape.v(3,22,k) = (-1)*1/4*(s^2-1)*(t-1)*t*(2*r+1);
    dshape.v(3,20,k) = (-1)*1/4*s*(s-1)*(t^2-1)*(2*r+1);
    dshape.v(3,24,k) = (-1)*1/4*(s^2-1)*t*(t+1)*(2*r+1);
    dshape.v(3,26,k) = (-1)*1/4*s*(s+1)*(t^2-1)*(2*r+1);
    dshape.v(3,5,k) = 1/2*(s^2-1)*(t^2-1)*(2*r-1);
    dshape.v(3,13,k) = (s^2-1)*(t-1)*t*r;
    dshape.v(3,17,k) = s*(s+1)*(t^2-1)*r;
    dshape.v(3,15,k) = (s^2-1)*t*(t+1)*r;
    dshape.v(3,11,k) = (s-1)*s*(t^2-1)*r;
    dshape.v(3,23,k) = 1/2*(s^2-1)*(t^2-1)*(2*r+1);
    dshape.v(3,14,k) = (-2)*(s^2-1)*(t^2-1)*r;
    
    % pressure shape functions
    shape.p(1,k) = 1/8*(1-s)*(1-t)*(1-r);
    shape.p(3,k) = 1/8*(1+s)*(1-t)*(1-r);
    shape.p(4,k) = 1/8*(1+s)*(1+t)*(1-r);
    shape.p(2,k) = 1/8*(1-s)*(1+t)*(1-r);
    shape.p(5,k) = 1/8*(1-s)*(1-t)*(1+r);
    shape.p(7,k) = 1/8*(1+s)*(1-t)*(1+r);
    shape.p(8,k) = 1/8*(1+s)*(1+t)*(1+r);
    shape.p(6,k) = 1/8*(1-s)*(1+t)*(1+r);
    
    % derivatives of pressure shape functions w.r.t s
    dshape.p(1,1,k) = (-1)*1/8*(t-1)*(r-1);
    dshape.p(1,3,k) = 1/8*(t-1)*(r-1);
    dshape.p(1,4,k) = (-1)*1/8*(t+1)*(r-1);
    dshape.p(1,2,k) = 1/8*(t+1)*(r-1);
    dshape.p(1,5,k) = 1/8*(t-1)*(r+1);
    dshape.p(1,7,k) = (-1)*1/8*(t-1)*(r+1);
    dshape.p(1,8,k) = 1/8*(t+1)*(r+1);
    dshape.p(1,6,k) = (-1)*1/8*(t+1)*(r+1);
    
    % derivatives of pressure shape functions w.r.t t
    dshape.p(2,1,k) = (-1)*1/8*(s-1)*(r-1);
    dshape.p(2,3,k) = 1/8*(s+1)*(r-1);
    dshape.p(2,4,k) = (-1)*1/8*(s+1)*(r-1);
    dshape.p(2,2,k) = 1/8*(s-1)*(r-1);
    dshape.p(2,5,k) = 1/8*(s-1)*(r+1);
    dshape.p(2,7,k) = (-1)*1/8*(s+1)*(r+1);
    dshape.p(2,8,k) = 1/8*(s+1)*(r+1);
    dshape.p(2,6,k) = (-1)*1/8*(s-1)*(r+1);
    
    % derivatives of pressure shape functions w.r.t r
    dshape.p(3,1,k) = (-1)*1/8*(s-1)*(t-1);
    dshape.p(3,3,k) = 1/8*(s+1)*(t-1);
    dshape.p(3,4,k) = (-1)*1/8*(s+1)*(t+1);
    dshape.p(3,2,k) = 1/8*(s-1)*(t+1);
    dshape.p(3,5,k) = 1/8*(s-1)*(t-1);
    dshape.p(3,7,k) = (-1)*1/8*(s+1)*(t-1);
    dshape.p(3,8,k) = 1/8*(s+1)*(t+1);
    dshape.p(3,6,k) = (-1)*1/8*(s-1)*(t+1);
end

end