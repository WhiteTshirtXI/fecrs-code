function [Vf,Vb] = rates_blm_fecrs(Y,epsilon,boundary_data)
% This code is for the case of planar deformation. Simple alterations can
% be made for a 3d scenario.
% blm: bead and link model.
% Boundary_data contains information about the boundary for calculation of
% the stokeslet interaction with the boundary.
% Boundary data should be given as a 3xN matrix (or else requires
% conversion).

%% reshape input data as required for bead model functions

N       = length(Y)/3;      % number of filament nodes
X(1,:)  = Y(1:N);
X(2,:)  = Y(N+1:2*N);
X(3,:)  = Y(2*N+1:end);
Y       = X;

M           = size(boundary_data,1)/6;  % number of boundary nodes on each face of cube.
b_data(1,:) = boundary_data(:,1);
b_data(2,:) = boundary_data(:,2);
b_data(3,:) = boundary_data(:,3);

%% main 

S   = 5e4;                  % non-dimensional parameter

N   = size(Y,2);            % number of beads
if size(Y,1) == 2
    Y = [Y; zeros(1,N)];    % augment Y with z component data if required
end

Q   = N-1;                  % number of links
ds  = 1/Q;                  % bead separation/equilibruim distance
    
F   = zeros(3,N);
Vf  = zeros(3,N);

F   = F + Q*get_bending_forces(Y) + S*get_spring_forces(Y,ds);
F   = reshape(F',[],1);

stokeslets_f    = get_reg_stokeslets(Y,Y,epsilon);                  % collocation with filament
Vf              = (1/8/pi)*stokeslets_f*F;                          % velocities at filament

stokeslets_b    = get_reg_stokeslets(Y,b_data,epsilon);             % collocation with boundary
Vb              = (1/8/pi)*stokeslets_b*F;                          % velocities at boundary due to filament

end % function


%% END OF FUNCTION
%  AUXILLARY FUNCTIONS BELOW

%% REGULARIZED STOKESLETS

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

%% BENDING FORCES

function F = get_bending_forces(x)
% Calculates restorative bending forces at Np 3-dimensional beads to keep a chain aligned

Np = size(x,2);
F = zeros(3,Np);

% bead 1.
bj   = x(:,1) - x(:,2);
bjp1 = x(:,2) - x(:,3);
F(:,1) = bjp1/norm(bj)/norm(bjp1) + dot(bj,bjp1)*(-bj)/norm(bj)^3/norm(bjp1);

% bead 2.
bjn1 = x(:,1) - x(:,2);
bj   = x(:,2) - x(:,3);
bjp1 = x(:,3) - x(:,4);
Fb = (bjn1 - bj)/norm(bjn1)/norm(bj) + dot(bjn1,bj)*(bjn1/norm(bjn1)^3/norm(bj) - bj/norm(bjn1)/norm(bj)^3);
Fc = bjp1/norm(bj)/norm(bjp1) + dot(bj,bjp1)*(-bj)/norm(bj)^3/norm(bjp1);
F(:,2) = Fb+Fc;

% beads 3 to Np-2.
for kk = 3:(Np-2)
    bjn2 = x(:,kk-2) - x(:,kk-1);
    bjn1 = x(:,kk-1) - x(:,kk);
    bj   = x(:,kk)   - x(:,kk+1);
    bjp1 = x(:,kk+1) - x(:,kk+2);
    Fa = -bjn2/norm(bjn2)/norm(bjn1)   + dot(bjn2,bjn1)*bjn1/norm(bjn2)/norm(bjn1)^3;
    Fb = (bjn1-bj)/norm(bjn1)/norm(bj) + dot(bjn1,bj)*(bjn1/norm(bjn1)^3/norm(bj) - bj/norm(bjn1)/norm(bj)^3);
    Fc = bjp1/norm(bj)/norm(bjp1)      + dot(bj,bjp1)*(-bj)/norm(bj)^3/norm(bjp1);
    F(:,kk) = Fa+Fb+Fc;
end

% bead Np-1.
bjn2 = x(:,Np-3) - x(:,Np-2);
bjn1 = x(:,Np-2) - x(:,Np-1);
bj   = x(:,Np-1) - x(:,Np);
Fa = -bjn2/norm(bjn2)/norm(bjn1)   + dot(bjn2,bjn1)*bjn1/norm(bjn2)/norm(bjn1)^3;
Fb = (bjn1-bj)/norm(bjn1)/norm(bj) + dot(bjn1,bj)*(bjn1/norm(bjn1)^3/norm(bj) - bj/norm(bjn1)/norm(bj)^3);
F(:,Np-1) = Fa+Fb;

% bead Np.
bjn2 = x(:,Np-2) - x(:,Np-1);
bjn1 = x(:,Np-1) - x(:,Np);
F(:,Np) = -bjn2/norm(bjn2)/norm(bjn1) + dot(bjn2,bjn1)*bjn1/norm(bjn2)/norm(bjn1)^3;

end % function

%% SPRING FORCES

function F = get_spring_forces(x, b0)
% Calculates the spring forces on a connected chain of beads in 3-dimensional space.

Np = size(x,2);
F  = zeros(3,Np);

for p = 1:Np
    if p > 1
        bpn1  = x(:,p-1) - x(:,p);
        Nbpn1 = norm(bpn1);
        F1    = -2*(Nbpn1 - b0)*(-bpn1/Nbpn1);
    else
        F1 = 0;
    end
    if p < Np
        bp   = x(:,p) - x(:,p+1);
        Nbp  = norm(bp);
        F2   = -2*(Nbp - b0)*(bp/Nbp);
    else
        F2 = 0;
    end
    F(:,p) = F1 + F2;
end
end % function
