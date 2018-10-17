% Calculates the fem velocity solution at any point in the domain
% Lx x Ly x Lz using shape functions
%
% [u_eval,v_eval,w_eval] = fem_evalVel(mesh,x,y,z,u,v,w)
%
% ========================================================================
% Inputs
% ========================================================================
% mesh = structure containing mesh parameters
% (pnt_x,pnt_y,pnt_z) = point(s) at which to evaluate solution
% (u,v,w) = velocity solutions from fem
%
% ========================================================================
% Ouputs
% ========================================================================
% (u_eval,v_eval,w_eval) = velocity solution at point(s) (x,y,z)

function [u_eval,v_eval,w_eval] = fem_evalVel(mesh,pnt_x,pnt_y,pnt_z,u,v,w)

% unpack mesh struct
x_dom = mesh.Lx;
y_dom = mesh.Ly;
z_dom = mesh.Lz;
Nx = mesh.Nx;
Ny = mesh.Ny;
Nz = mesh.Nz;

nNodesDir = Nx+(Nx-1); % number of velociy grid point per direction

% sort parameters
pnt_x = pnt_x(:); pnt_y = pnt_y(:); pnt_z = pnt_z(:);

%% griddedInterpolant method

x = linspace(x_dom(1),x_dom(2),Nx+(Nx-1));
y = linspace(y_dom(1),y_dom(2),Ny+(Ny-1));
z = linspace(z_dom(1),z_dom(2),Nz+(Nz-1));

[x,y,z] = ndgrid(x,y,z);

u = reshape(u,[nNodesDir,nNodesDir,nNodesDir]);
u = permute(u, [2,1,3]);

v = reshape(v,[nNodesDir,nNodesDir,nNodesDir]);
v = permute(v, [2,1,3]);

w = reshape(w,[nNodesDir,nNodesDir,nNodesDir]);
w = permute(w, [2,1,3]);

F_u = griddedInterpolant(x,y,z,u);
F_v = griddedInterpolant(x,y,z,v);
F_w = griddedInterpolant(x,y,z,w);

for i=1:length(pnt_x)
    u_eval(i,1) = F_u(pnt_x(i),pnt_y(i),pnt_z(i));
    v_eval(i,1) = F_v(pnt_x(i),pnt_y(i),pnt_z(i));
    w_eval(i,1) = F_w(pnt_x(i),pnt_y(i),pnt_z(i));
end

end