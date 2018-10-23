% Constructs load vector and applies boundary conditions
%
% [b] = fem_calcLoadVec(mesh,u_bound,v_bound,w_bound)
%
% ========================================================================
% Inputs
% ========================================================================
% mesh = structure containing mesh parameters
% u_bound = structure containing boundary conditions on u
% v_bound = structure containing boundary conditions on v
% w_bound = structure containing boundary conditions on w
%
% ========================================================================
% Outputs
% ========================================================================
% b = load vector

function [b] = fem_calcLoadVec(mesh,Vb)

npf = length(Vb)/3/6;

% unpack mesh struct
Lx = mesh.Lx;
Ly = mesh.Ly;
Lz = mesh.Lz;
vCoords = mesh.vCoords;
nvN = mesh.nvN;

% unpack boundary condition structs
u_bound_face_1 = -Vb(1:npf);
u_bound_face_2 = -Vb(npf+1:2*npf);
u_bound_face_3 = -Vb(2*npf+1:3*npf);
u_bound_face_4 = -Vb(3*npf+1:4*npf);
u_bound_face_5 = -Vb(4*npf+1:5*npf);
u_bound_face_6 = -Vb(5*npf+1:6*npf);

v_bound_face_1 = -Vb(6*npf+1:7*npf);
v_bound_face_2 = -Vb(7*npf+1:8*npf);
v_bound_face_3 = -Vb(8*npf+1:9*npf);
v_bound_face_4 = -Vb(9*npf+1:10*npf);
v_bound_face_5 = -Vb(10*npf+1:11*npf);
v_bound_face_6 = -Vb(11*npf+1:12*npf);

w_bound_face_1 = -Vb(12*npf+1:13*npf);
w_bound_face_2 = -Vb(13*npf+1:14*npf);
w_bound_face_3 = -Vb(14*npf+1:15*npf);
w_bound_face_4 = -Vb(15*npf+1:16*npf);
w_bound_face_5 = -Vb(16*npf+1:17*npf);
w_bound_face_6 = -Vb(17*npf+1:18*npf);

% create load vector
b = zeros(3*mesh.nvN+mesh.npN,1);

% apply boundary conditions to load vector
ind_1 = find(vCoords(:,2) == Ly(2)); % indices for back face
ind_2 = find(vCoords(:,1) == Lx(2)); % indices for right face
ind_3 = find(vCoords(:,2) == Ly(1)); % indices for front face
ind_4 = find(vCoords(:,1) == Lx(1)); % indices for left face
ind_5 = find(vCoords(:,3) == Lz(1)); % indices for bottom face
ind_6 = find(vCoords(:,3) == Lz(2)); % indices for top face

b_u = b(1:nvN); % load vector corresponding to u
b_v = b(1+nvN:2*nvN); % load vector corresponding to v
b_w = b(2*nvN+1:3*nvN); % load vector corresponding to w
b_p = b((3*nvN+1):end); % load vector corresponding to p

b_u(ind_1) = u_bound_face_1;
b_v(ind_1) = v_bound_face_1;
b_w(ind_1) = w_bound_face_1;

b_u(ind_2) = u_bound_face_2;
b_v(ind_2) = v_bound_face_2;
b_w(ind_2) = w_bound_face_2;

b_u(ind_3) = u_bound_face_3;
b_v(ind_3) = v_bound_face_3;
b_w(ind_3) = w_bound_face_3;

b_u(ind_4) = u_bound_face_4;
b_v(ind_4) = v_bound_face_4;
b_w(ind_4) = w_bound_face_4;

b_u(ind_5) = u_bound_face_5;
b_v(ind_5) = v_bound_face_5;
b_w(ind_5) = w_bound_face_5;

b_u(ind_6) = u_bound_face_6;
b_v(ind_6) = v_bound_face_6;
b_w(ind_6) = w_bound_face_6;

b = [b_u; b_v; b_w; b_p]; % reform load vector

b = [b; 0]; % fix pressure

end
