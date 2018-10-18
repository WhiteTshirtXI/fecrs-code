% Constructs stiffness matrix and applies boundary conditions
%
% [M,mesh] = fem_calcStiffnessMat(Nx,Ny,Nz,Lx,Ly,Lz)
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
% M = stiffness matrix
% mesh = structure containing mesh parameters

function [M,mesh] = fem_calcStiffnessMat(Nx,Ny,Nz,Lx,Ly,Lz)

[GQ] = fem_setupQuadrature();
[shape, dshape] = fem_calcShapeFun(GQ);
[mesh] = fem_genMesh(Nx,Ny,Nz,Lx,Ly,Lz);
[elem] = fem_calcJacob(mesh,dshape);
[M] = fem_assemble(mesh,elem,GQ,shape);

% apply boundary conditions to stiffness matrix
vCoords = mesh.vCoords;
nvN = mesh.nvN;

ind_1 = find(vCoords(:,2) == Ly(2)); % indices for back face
ind_2 = find(vCoords(:,1) == Lx(2)); % indices for right face
ind_3 = find(vCoords(:,2) == Ly(1)); % indices for front face
ind_4 = find(vCoords(:,1) == Lx(1)); % indices for left face
ind_5 = find(vCoords(:,3) == Lz(1)); % indices for bottom face
ind_6 = find(vCoords(:,3) == Lz(2)); % indices for top face

M_u = M(1:nvN,:); % stiffness matrix corresponding to u
M_v = M(1+nvN:2*nvN,:); % stiffness matrix corresponding to v
M_w = M(2*nvN+1:3*nvN,:); % stiffness matrix corresponding to w
M_p = M((3*nvN+1):end,:); % stiffness matrix corresponding to p

M_u([ind_1, ind_2, ind_3, ind_4, ind_5, ind_6],:) = 0;
M_v([ind_1, ind_2, ind_3, ind_4, ind_5, ind_6],:) = 0;
M_w([ind_1, ind_2, ind_3, ind_4, ind_5, ind_6],:) = 0;

ind1 = sub2ind(size(M_u), ind_1, ind_1);
ind2 = sub2ind(size(M_v), ind_1, nvN+ind_1);
ind3 = sub2ind(size(M_w), ind_1, 2*nvN+ind_1);
ind4 = sub2ind(size(M_u), ind_2, ind_2);
ind5 = sub2ind(size(M_v), ind_2, nvN+ind_2);
ind6 = sub2ind(size(M_w), ind_2, 2*nvN+ind_2);
ind7 = sub2ind(size(M_u), ind_3, ind_3);
ind8 = sub2ind(size(M_v), ind_3, nvN+ind_3);
ind9 = sub2ind(size(M_w), ind_3, 2*nvN+ind_3);
ind10 = sub2ind(size(M_u), ind_4, ind_4);
ind11 = sub2ind(size(M_v), ind_4, nvN+ind_4);
ind12 = sub2ind(size(M_w), ind_4, 2*nvN+ind_4);
ind13 = sub2ind(size(M_u), ind_5, ind_5);
ind14 = sub2ind(size(M_v), ind_5, nvN+ind_5);
ind15 = sub2ind(size(M_w), ind_5, 2*nvN+ind_5);
ind16 = sub2ind(size(M_u), ind_6, ind_6);
ind17 = sub2ind(size(M_v), ind_6, nvN+ind_6);
ind18 = sub2ind(size(M_w), ind_6, 2*nvN+ind_6);

M_u([ind1,ind4,ind7,ind10,ind13,ind16]) = 1;
M_v([ind2,ind5,ind8,ind11,ind14,ind17]) = 1;
M_w([ind3,ind6,ind9,ind12,ind15,ind18]) = 1;

M = [M_u; M_v; M_w; M_p]; % reform stiffness matrix

M = [M; zeros(1,length(M))]; M(end,end) = 1; % pin pressure

end