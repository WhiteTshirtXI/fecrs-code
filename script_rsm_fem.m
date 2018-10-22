% This script uses the elastohydrodynamic model using the
% method of regularized stokeslets along with the finite element
% method (for Stokes flow), to simulate a relaxing filament in a
% bounded rectangular domain.

clear all; close all;
tic;

%% options
save_data           = 0;                    % choose whether to save data
filename            = 'data_blm_fem.mat';   % filename for save
nSave               = 10000;                % save every nSave timesteps

%% fem parameters

% nodes per direction
nx                  = 5;
ny                  = 5;
nz                  = 5;

% domain size
lx                  = [0,1];
ly                  = [0,1];
lz                  = [0,1];

%% other parameters

% elastohydrodynamic model w/ reg. stokeslets
eps                 = 0.01;                % epsilon
a                   = 0.5;                 % for the initial config y=ax^2
Q                   = 20;                  % number of segments
N                   = Q+1;                 % number of joints/beads

% time setup
tStart              = 0;
tEnd                = 1e-3;
tSpan               = linspace(tStart,tEnd,25);
Nt                  = length(tSpan);
plotcolors          = [linspace(0.5,0,Nt)',linspace(1,0,Nt)',ones(Nt,1)];

%% setup

% generate initial configuration
[x, y, ~] = get_parabolic_nodes(a,Q);

% transform bead coords to domain [0,1]^3 - need to generalise
x = x + 0.5;
y = y + 0.25;

X = [x',y',0.5*ones(N,1)];      % add z component
s = calc_arclength(X(:,1:2)');  % calculates arclength (check)

%% find boundary nodes

% number of velocity nodes per direction
nxv = nx+(nx-1);
nyv = ny+(ny-1);
nzv = nz+(nz-1);

% velocity node coordinates
x = linspace(lx(1),lx(2),nxv);
y = linspace(ly(1),ly(2),nyv);
z = linspace(lz(1),lz(2),nzv);

[x,y,z] = meshgrid(x,y,z);
x = x(:); y = y(:); z = z(:); % node coordinates

% boundary coordinate indices
ind1 = find(y == ly(2));
ind2 = find(x == lx(2));
ind3 = find(y == ly(1));
ind4 = find(x == lx(1));
ind5 = find(z == lz(1));
ind6 = find(z == lz(2));

% coordinates of boundary nodes
bnd = [x(ind1),y(ind1),z(ind1);x(ind2),y(ind2),z(ind2); ...
       x(ind3),y(ind3),z(ind3);x(ind4),y(ind4),z(ind4); ...
       x(ind5),y(ind5),z(ind5);x(ind6),y(ind6),z(ind6)];

%% solve problem

% calculate stiffness matrix for fem
[M,mesh] = fem_calcStiffnessMat(nx,ny,nz,lx,ly,lz);

% calculate inverse
M       = full(M);
Minv    = pinv(M);

Y0      = [X(:,1);X(:,2);X(:,3)]; 
tspan   = linspace(tStart,tEnd,30);

dY_fem_rsm      = @(t,Y) rates_fem_rsm(t,Y,Minv,mesh,bnd,eps,save_data,nSave);
sol_struct      = ode15s(dY_fem_rsm,[tStart,tEnd],Y0);
sol_fem_rsm     = deval(sol_struct,tSpan);
sol_dt          = diff(sol_struct.x);
sol_dt_sum      = cumsum(sol_dt);

fprintf('   Solve complete!\n')

%% plot initial and final configs
figure(1);
for kk = 1:Nt
    subplot(1,2,1); hold on; box on;
    x_end = sol_fem_rsm(1:N,kk);
    y_end = sol_fem_rsm(N+1:2*N,kk);
    h(kk) = plot(x_end,y_end,'.-','Color',plotcolors(kk,:));
    xlim([-0.1,1.1])
    ylim([0.1,0.6])
    axis equal
end
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
legend([h(1), h(end)],{'int. config','final config'})

figure(1);
subplot(1,2,2); box on;
plot(sol_dt_sum,sol_dt)
xlabel('$t$','Interpreter','latex')
ylabel('$\delta t$','Interpreter','latex')

%% completion
script_walltime = toc;
fprintf('   Script complete in %g seconds!\n',script_walltime)

