% This scipt uses the bead and link model along with the finite element
% method (for Stokes flow), to simulate a relaxing filament in a
% bounded rectangular domain.

clear all;

%% options
save_data           = 0;                   % choose whether to save data
filename            = 'data_blm_fem.mat';  % filename for save
nSave               = 10000;                 % save every nSave timesteps

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

% bead and link model
eps                 = 0.01;                % epsilon
a                   = 0.5;                 % for the initial config y=ax^2
Q                   = 60;                  % number of segments

% time stepping
tStart              = 0;
tEnd                = 1e-2;
dt                  = 1e-7;

%% setup

% time stepping
tvec = tStart:dt:tEnd;
Nt = length(tvec);

% generate initial configuration
[x, y, ~] = get_parabolic_nodes(a,Q);

% transform bead coords to domain [0,1]^3 - need to generalise
x = x + 0.5;
y = y + 0.25;

nBeads = length(x);

X = [x',y',0.5*ones(nBeads,1)]; % add z component

s = calc_arclength(X(:,1:2)'); % calculates arclength

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

count = 1;

% calculate stiffness matrix for fem
[M,mesh] = fem_calcStiffnessMat(nx,ny,nz,lx,ly,lz);

% calculate inverse
M = full(M);
Minv = pinv(M);

for t=1:Nt
    
    beads = X(:,:,t);
    
    % solve blm for bead and boundary velocities
    [Vf,Vb] = rates_blm_fecrs(beads(:),eps,bnd);
    
    % calculate load vector for fem
    [b] = fem_calcLoadVec(mesh,Vb);
    
    us = Vf(1:(Q+1));
    vs = Vf(Q+2:(2*Q+2));
    ws = Vf(2*Q+3:(3*Q+3));
    
    % solve system
    c = Minv*b;
    
    % velocity solutions at nodes of fem mesh
    U(:,t) = c(1:mesh.nvN); % u solution
    V(:,t) = c(1+mesh.nvN:2*mesh.nvN); % v solution
    W(:,t) = c(2*mesh.nvN+1:3*mesh.nvN); % w solution
    %P(:,t) = c((3*mesh.nvN+1):end); % p solution
    
    % evaluate fem solution at all beads
    [U_bead,V_bead,W_bead] = fem_evalVel(mesh,X(:,1,t),X(:,2,t),X(:,3,t)...
        ,U(:,t),V(:,t),W(:,t));
    
    % evaluate final velocity at beads (u = U + us)
    for i=1:(Q+1)
        u_bead(i,t) = U_bead(i) + us(i);
        v_bead(i,t) = V_bead(i) + vs(i);
        w_bead(i,t) = W_bead(i) + ws(i);
    end
    
    % time step - forward Euler
    for i=1:(Q+1)
        X(i,:,t+1) = timestep_forwardEuler(u_bead(i,t),v_bead(i,t),...
            w_bead(i,t),X(i,:,t),dt);
    end
    
    s(t+1) = calc_arclength(X(:,1:2,t)'); % calculate arclength
    
    perccount(t,Nt); % displays percentage complete
    
   
    if (save_data == 1) 
        if (mod(t,nSave) == 0) % save every nSave timesteps
            Xsave(:,:,count) = X(:,:,t);
            save(filename,'Xsave','tvec','nBeads','u_bead','v_bead','w_bead','s'); 
            count = count + 1;
        end
    end
    
end

%% plot initial and final configs
figure;
scatter(X(:,1,1),X(:,2,1),'.r');
hold on;
scatter(X(:,1,end),X(:,2,end),'.b');
xlim([0,1])
ylim([0,1])
axis equal
