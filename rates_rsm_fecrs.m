function [Vf,Vb]=rates_rsm_fecrs(Y,epsilon,boundary_data)

% this code (previously called CalcRates.m) calculates the movement of a
% filament using the EHD model with forces decribed by the method of
% regularised stokeslets. Code by Smith (2018).

% optional methods of calculating blocks have been removed for brevity of
% the code. Please check "rates_rsm.m" to view these omissions.

% this code is formatted such that output data Vf and Vb contain the
% filament and boundary velocities respectively, to be fed into the finite
% element method (FEM) by C.V. Neal (2017).

%% check dimensions of Y and transpose accordingly
%  ode15s takes 1xQ input Y and transposes to Qx1
%  my manual solver does not do this, so we need to check and transpose

if size(Y,1)==1
    Y = Y';
end

%% calculate thetas using interpolating spline
%  x,y,z,th should be row vectors!

N   = size(Y,1)/3;
x   = Y(1:N)';
y   = Y(N+1:2*N)';
z   = Y(2*N+1:3*N)';
X   = [x;y;z];
bX  = diff(X')';

Q       = N-1;
ds      = 1/Q;
DS      = sqrt(diff(x).^2 + diff(y).^2);
s       = [0,cumsum(DS)];
delta   = max(s)/Q;
ss      = 0:delta:max(s);

yy      = spline(s,y,ss);
xx      = spline(s,x,ss);
dy      = diff(yy);
dx      = diff(xx);
th      = atan(dy./dx);

%% setup transformation matrix for hydrodynamics calculation

R = [ cos(th)     -sin(th)       zeros(1,Q)   ;
      sin(th)      cos(th)       zeros(1,Q)   ;
      zeros(1,Q)   zeros(1,Q)    ones(1,Q)     ] ;

%% segment midpoints
xm = (x(1:end-1)+x(2:end))/2;
ym = (y(1:end-1)+y(2:end))/2;
zm = (z(1:end-1)+z(2:end))/2;
Xm      = [xm;ym;zm];

%% hydrodynamics block: velocity-force

AH3 = -RegStokesletAnalyticIntegrals(Xm',Xm',ds/2,R,epsilon);
AH2 = AH3(1:2*Q,1:2*Q); % just extract (x,y) components for planar motion problem

%% kinematic block:

AK = [ones(Q,1)  zeros(Q,1) tril(repmat(-ds*sin(th(:))',Q,1),-1)+diag(-ds/2*sin(th(:))') ;
      zeros(Q,1)  ones(Q,1) tril(repmat( ds*cos(th(:))',Q,1),-1)+diag( ds/2*cos(th(:))')   ];

%% elastodynamics block: moment-force

AE=-[triu(repmat( ds*y(1:Q)',1,Q))+ triu(repmat(-ds*ym,Q,1)) ...
     triu(repmat(-ds*x(1:Q)',1,Q))+ triu(repmat( ds*xm,Q,1)) ;
     ds*ones(1,Q)                   zeros(1,Q)               ;
     zeros(1,Q)                     ds*ones(1,Q)            ];

%% construct linear system and solve

A       = [zeros(Q+2) AE; AK AH2]; 
dtheta  = diff(th)'/ds; 

b = [  0                      ;      % zero moment about x{1}
      dtheta                  ;      % elastohydrodynamic moment balance rows about x{2}...x{N+1}
      zeros(2,1)              ;      % total force balance rows
      zeros(2*Q,1)              ]  ; % velocity rows
  
dZ      = A\b;
forces  = [dZ(Q+3:end);zeros(Q,1)];  % zero force in z due to planar motion.

%% calculate filament bead/joint velocities

stokeslets_f    = get_reg_stokeslets(Xm,X,epsilon);
Vf              = (1/8/pi)*stokeslets_f*forces; 

%% calculate boundary velocities

M           = size(boundary_data,1)/6;  % number of boundary nodes on each face of cube.
b_data(1,:) = boundary_data(:,1);
b_data(2,:) = boundary_data(:,2);
b_data(3,:) = boundary_data(:,3);

stokeslets_b    = get_reg_stokeslets(Xm,b_data,epsilon);
Vb              = (1/8/pi)*stokeslets_b*forces;

%% test

figure(10); subplot(1,2,1); 
cla
quiver(xm,ym,forces(1:Q)',forces(Q+1:2*Q)')
hold on
quiver(x,y,Vf(1:N)',Vf(N+1:2*N)')
plot(x,y,'o-')
axis equal
xlabel('x')
ylabel('y')

figure(10); subplot(1,2,2); cla
plot(th)

drawnow


end % function

%% END OF MAIN FUNCTION
%  AUXILLARY FUNCTIONS ARE LISTED BELOW

%% ANALYTIC REG STOKESLETS

function D=RegStokesletAnalyticIntegrals(x,X,h,R,epsilon)

% x is 3Nx1 field points
% X is 3Qx1 centres of intervals of integration
% h is half-length of integration
% R is a 3x3Q rotation matrix into the local coordinates of each interval
% of integration

x=x(:);
X=X(:);
[x1,x2,x3]=ExtractComponents(x);
[X1,X2,X3]=ExtractComponents(X);

N=length(x1);
Q=length(X1);

RM=kron(R,ones(N,1)); % RM is 3N x 3Q. Repeat for every field point.

RM11=repmat(R(1,      1:Q),N,1);
RM12=repmat(R(1,  Q+1:2*Q),N,1);
RM13=repmat(R(1,2*Q+1:3*Q),N,1);

RM21=repmat(R(2,      1:Q),N,1);
RM22=repmat(R(2,  Q+1:2*Q),N,1);
RM23=repmat(R(2,2*Q+1:3*Q),N,1);

RM31=repmat(R(3,      1:Q),N,1);
RM32=repmat(R(3,  Q+1:2*Q),N,1);
RM33=repmat(R(3,2*Q+1:3*Q),N,1);

xi1=x1*ones(1,Q)-ones(N,1)*X1';
xi2=x2*ones(1,Q)-ones(N,1)*X2';
xi3=x3*ones(1,Q)-ones(N,1)*X3';

% rotate into local coordinates.
xiL1=RM11.*xi1 + RM21.*xi2 + RM31.*xi3;
xiL2=RM12.*xi1 + RM22.*xi2 + RM32.*xi3;
xiL3=RM13.*xi1 + RM23.*xi2 + RM33.*xi3;

denom=xiL2.^2+xiL3.^2+epsilon^2;
idenom=1./denom;
repsM=sqrt((xiL1-h).^2+denom);
repsP=sqrt((xiL1+h).^2+denom);
irepsM=1./repsM;
irepsP=1./repsP;

I11=(-(xiL1-h).*irepsM+(xiL1+h).*irepsP).*(epsilon^2.*idenom-1)...
    +2*log((h-xiL1+repsM)./(-h-xiL1+repsP));
I22=(-(xiL1-h).*irepsM+(xiL1+h).*irepsP).*(epsilon^2+xiL2.^2).*idenom...
      +log((h-xiL1+repsM)./(-h-xiL1+repsP));
I33=(-(xiL1-h).*irepsM+(xiL1+h).*irepsP).*(epsilon^2+xiL3.^2).*idenom...
      +log((h-xiL1+repsM)./(-h-xiL1+repsP));
I12=xiL2.*(irepsM-irepsP);
I21=I12;
I13=xiL3.*(irepsM-irepsP);
I31=I13;
I23=(-(xiL1-h).*irepsM+(xiL1+h).*irepsP).*xiL2.*xiL3.*idenom;
I32=I23;

% transform 2nd rank tensor - might be neater and safer with cell arrays
D11=(RM11.*I11+RM12.*I21+RM13.*I31).*RM11...
   +(RM11.*I12+RM12.*I22+RM13.*I32).*RM12...
   +(RM11.*I13+RM12.*I23+RM13.*I33).*RM13;

D12=(RM11.*I11+RM12.*I21+RM13.*I31).*RM21...
   +(RM11.*I12+RM12.*I22+RM13.*I32).*RM22...
   +(RM11.*I13+RM12.*I23+RM13.*I33).*RM23;

D13=(RM11.*I11+RM12.*I21+RM13.*I31).*RM31...
   +(RM11.*I12+RM12.*I22+RM13.*I32).*RM32...
   +(RM11.*I13+RM12.*I23+RM13.*I33).*RM33;

D21=(RM21.*I11+RM22.*I21+RM23.*I31).*RM11...
   +(RM21.*I12+RM22.*I22+RM23.*I32).*RM12...
   +(RM21.*I13+RM22.*I23+RM23.*I33).*RM13;

D22=(RM21.*I11+RM22.*I21+RM23.*I31).*RM21...
   +(RM21.*I12+RM22.*I22+RM23.*I32).*RM22...
   +(RM21.*I13+RM22.*I23+RM23.*I33).*RM23;

D23=(RM21.*I11+RM22.*I21+RM23.*I31).*RM31...
   +(RM21.*I12+RM22.*I22+RM23.*I32).*RM32...
   +(RM21.*I13+RM22.*I23+RM23.*I33).*RM33;

D31=(RM31.*I11+RM32.*I21+RM33.*I31).*RM11...
   +(RM31.*I12+RM32.*I22+RM33.*I32).*RM12...
   +(RM31.*I13+RM32.*I23+RM33.*I33).*RM13;

D32=(RM31.*I11+RM32.*I21+RM33.*I31).*RM21...
   +(RM31.*I12+RM32.*I22+RM33.*I32).*RM22...
   +(RM31.*I13+RM32.*I23+RM33.*I33).*RM23;

D33=(RM31.*I11+RM32.*I21+RM33.*I31).*RM31...
   +(RM31.*I12+RM32.*I22+RM33.*I32).*RM32...
   +(RM31.*I13+RM32.*I23+RM33.*I33).*RM33;

D=[D11 D12 D13; D21 D22 D23; D31 D32 D33]/8/pi;

end % function

%% EXTRACT PLANAR COMPONENTS

function [x1,x2,x3]=ExtractComponents(x)

  N=length(x)/3;
  x1=x(1:N);
  x2=x(N+1:2*N);
  x3=x(2*N+1:3*N);

end % function

%% GET REGULARIZED STOKESLETS

function S = get_reg_stokeslets(x, x0, epsilon)
% Calculates the Oseen tensor using the method of regularised stokeslets.
% x is the 3xM matrix of stokeslet positions
% x0 is the 3xN matrix of positions at which velocity due to the stokeslets
% at positions x is to be calculated.

eps2 = epsilon^2;

N   = size(x,2);
M   = size(x0,2);

Sxx = zeros(M,N); Sxy = zeros(M,N); Sxz = zeros(M,N);
Syx = zeros(M,N); Syy = zeros(M,N); Syz = zeros(M,N);
Szx = zeros(M,N); Szy = zeros(M,N); Szz = zeros(M,N);

for p = 1:M

    rx = x(1,:) - x0(1,p);
    ry = x(2,:) - x0(2,p);
    rz = x(3,:) - x0(3,p);

    r2 = rx.^2 + ry.^2 + rz.^2;
    r_eps = (r2 + eps2).^1.5;
    r_num = (r2 + 2*eps2);

    Sxx(p,:) = (r_num + rx.*rx)./r_eps;
    Sxy(p,:) = (rx.*ry)./r_eps;
    Sxz(p,:) = (rx.*rz)./r_eps;

    Syx(p,:) = (ry.*rx)./r_eps;
    Syy(p,:) = (r_num + ry.*ry)./r_eps;
    Syz(p,:) = (ry.*rz)./r_eps;

    Szx(p,:) = (rz.*rx)./r_eps;
    Szy(p,:) = (rz.*ry)./r_eps;
    Szz(p,:) = (r_num + rz.*rz)./r_eps;

end

A1 = horzcat(Sxx, Sxy, Sxz);
A2 = horzcat(Syx, Syy, Syz);
A3 = horzcat(Szx, Szy, Szz);
S  = vertcat(A1, A2, A3);

end % function
