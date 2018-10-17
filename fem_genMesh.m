% Generates a hexahedron mesh on the rectangular domain Lx x Ly x Lz
%
% [mesh] = fem_genMesh(Nx,Ny,Nz,Lx,Ly,Lz)
%
% ========================================================================
% Inputs
% ========================================================================
% Nx = number of grid points in the x direction
% Ny = number of grid points in the y direction
% Nz = number of grid points in the z direction
% Lx = domain of x
% Ly = domain of y
% Lz = domain of z
%
% ========================================================================
% Outputs
% ========================================================================
% mesh.nE = number of elements in mesh
% mesh.nvN = number of total velocity nodes in mesh
% mesh.npN = number of total pressure nodes in mesh
% mesh.vCoords  = coordinates of velocity nodes
% mesh.pCoords = coordinates of pressure nodes
% mesh.vConnect = connectivity matrix for velocity nodes
% mesh.pConnect = connectivity matrix for pressure nodes
% mesh.Nx = number of grid points in the x direction
% mesh.Ny = number of grid points in the y direction
% mesh.Nz = number of grid points in the z direction
% mesh.Lx = domain of x
% mesh.Ly = domain of y
% mesh.Lz = domain of z

function [mesh] = fem_genMesh(Nx,Ny,Nz,Lx,Ly,Lz)

mesh.Lx = Lx;
mesh.Ly = Ly;
mesh.Lz = Lz;
mesh.Nx = Nx;
mesh.Ny = Ny;
mesh.Nz = Nz;

% mesh properties
mesh.nE = (Nx-1)*(Ny-1)*(Nz-1); % number of elements
mesh.nvN = (Nx+(Nx-1))*(Ny+(Ny-1))*(Nz+(Nz-1)); % number of v nodes
mesh.npN = Nx*Ny*Nz; % number of p nodes

% nodes for pressure grid
xp = linspace(Lx(1),Lx(2),Nx);
yp = linspace(Ly(1),Ly(2),Ny);
zp = linspace(Lz(1),Lz(2),Nz);

[xp,yp,zp] = meshgrid(xp,yp,zp);
mesh.pCoords = [xp(:) yp(:) zp(:)];

% nodes for velocity grid
xv = linspace(Lx(1),Lx(2),Nx+(Nx-1));
yv = linspace(Ly(1),Ly(2),Ny+(Ny-1));
zv = linspace(Lz(1),Lz(2),Nz+(Nz-1));

[xv,yv,zv] = meshgrid(xv,yv,zv);
mesh.vCoords = [xv(:) yv(:) zv(:)];

% generate velocity connectivity matrix
mesh.vConnect = zeros(mesh.nE,27);
m = 1;
for j=1:2:(Ny+(Ny-1))
    for i=1:2*(Ny+(Ny-1)):(Nx+(Nx-1))*(Ny+(Ny-1))
        for k=1:2*(Nx+(Nx-1))*(Ny+(Ny-1)):mesh.nvN
            % check if node is on boundary Lx(2), Ly(2) or Lz(2)
            if (mesh.vCoords(i,1) == Lx(2) || mesh.vCoords(j,2) == ...
                    Ly(2) || mesh.vCoords(k,3) == Lz(2))
                % do not update matrix
            else
                % update matrix with nodal connectivity
                mesh.vConnect(m,1) = i+(j-1)+(k-1);
                mesh.vConnect(m,2) = i+1+(j-1)+(k-1);
                mesh.vConnect(m,3) = i+2+(j-1)+(k-1);
                mesh.vConnect(m,4) = (Ny+(Ny-1))+i+(j-1)+(k-1);
                mesh.vConnect(m,5) = (Ny+(Ny-1))+i+1+(j-1)+(k-1);
                mesh.vConnect(m,6) = (Ny+(Ny-1))+i+2+(j-1)+(k-1);
                mesh.vConnect(m,7) = 2*(Ny+(Ny-1))+i+(j-1)+(k-1);
                mesh.vConnect(m,8) = 2*(Ny+(Ny-1))+ i+1+(j-1)+(k-1);
                mesh.vConnect(m,9) = 2*(Ny+(Ny-1))+i+2+(j-1)+(k-1);
                mesh.vConnect(m,10) = (Nx+(Nx-1))*(Ny+(Ny-1))+i+...
                    (j-1)+(k-1);
                mesh.vConnect(m,11) = (Nx+(Nx-1))*(Ny+(Ny-1))+i+1+...
                    (j-1)+(k-1);
                mesh.vConnect(m,12) = (Nx+(Nx-1))*(Ny+(Ny-1))+i+2+...
                    (j-1)+(k-1);
                mesh.vConnect(m,13) = (Nx+(Nx-1))*(Ny+(Ny-1))+(Ny+...
                    (Ny-1))+i+(j-1)+(k-1);
                mesh.vConnect(m,14) = (Nx+(Nx-1))*(Ny+(Ny-1))+(Ny+...
                    (Ny-1))+i+1+(j-1)+(k-1);
                mesh.vConnect(m,15) = (Nx+(Nx-1))*(Ny+(Ny-1))+(Ny+...
                    (Ny-1))+i+2+(j-1)+(k-1);
                mesh.vConnect(m,16) = (Nx+(Nx-1))*(Ny+(Ny-1))+2*(Ny+...
                    (Ny-1))+i+(j-1)+(k-1);
                mesh.vConnect(m,17) = (Nx+(Nx-1))*(Ny+(Ny-1))+2*(Ny+...
                    (Ny-1))+i+1+(j-1)+(k-1);
                mesh.vConnect(m,18) = (Nx+(Nx-1))*(Ny+(Ny-1))+2*(Ny+...
                    (Ny-1))+i+2+(j-1)+(k-1);
                mesh.vConnect(m,19) = 2*(Nx+(Nx-1))*(Ny+(Ny-1))+i+...
                    (j-1)+(k-1);
                mesh.vConnect(m,20) = 2*(Nx+(Nx-1))*(Ny+(Ny-1))+i+1+...
                    (j-1)+(k-1);
                mesh.vConnect(m,21) = 2*(Nx+(Nx-1))*(Ny+(Ny-1))+i+2+...
                    (j-1)+(k-1);
                mesh.vConnect(m,22) = 2*(Nx+(Nx-1))*(Ny+(Ny-1))+(Ny+...
                    (Ny-1))+i+(j-1)+(k-1);
                mesh.vConnect(m,23) = 2*(Nx+(Nx-1))*(Ny+(Ny-1))+(Ny+...
                    (Ny-1))+i+1+(j-1)+(k-1);
                mesh.vConnect(m,24) = 2*(Nx+(Nx-1))*(Ny+(Ny-1))+(Ny+...
                    (Ny-1))+i+2+(j-1)+(k-1);
                mesh.vConnect(m,25) = 2*(Nx+(Nx-1))*(Ny+(Ny-1))+2*(Ny+...
                    (Ny-1))+i+(j-1)+(k-1);
                mesh.vConnect(m,26) = 2*(Nx+(Nx-1))*(Ny+(Ny-1))+2*(Ny+...
                    (Ny-1))+i+1+(j-1)+(k-1);
                mesh.vConnect(m,27) = 2*(Nx+(Nx-1))*(Ny+(Ny-1))+2*(Ny+...
                    (Ny-1))+i+2+(j-1)+(k-1);
                m = m+1;
            end
        end
    end
end

% generate pressure connectivity matrix
mesh.pConnect = zeros(mesh.nE,8);
m = 1;
for j=1:Ny
    for i=1:Ny:Nx*Ny
        for k=1:(Nx*Ny):mesh.npN
            % check if node is on boundary Lx(2), Ly(2) or Lz(2)
            if (mesh.pCoords(i,1) == Lx(2) || mesh.pCoords(j,2) == Ly(2)...
                    || mesh.pCoords(k,3) == Lz(2))
                % do not update matrix
            else
                mesh.pConnect(m,1) = i+(j-1)+(k-1);
                mesh.pConnect(m,2) = i+1+(j-1)+(k-1);
                mesh.pConnect(m,3) = Ny+i+(j-1)+(k-1);
                mesh.pConnect(m,4) = Ny+i+1+(j-1)+(k-1);
                mesh.pConnect(m,5) = Nz*Ny+i+(j-1)+(k-1);
                mesh.pConnect(m,6) = Nz*Ny+i+1+(j-1)+(k-1);
                mesh.pConnect(m,7) = Nz*Ny+Ny+ i+(j-1)+(k-1);
                mesh.pConnect(m,8) = Nz*Ny+Ny+i+1+(j-1)+(k-1);
                m=m+1;
            end
        end
    end
end