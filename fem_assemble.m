% Assembles full stiffness matrix using element matrices
%
% [M] = fem_assemble(mesh,elem,GQ,shape)
%
% ========================================================================
% Inputs
% ========================================================================
% mesh = structure containing mesh parameters
% elem = structure consisting of jacobian matrices, element matrices e.t.c
% GQ = structure containing required Gaussian quadrature nodes and weights
% shape = structure containing velocity and pressure shape functions
% evaluated at GQ nodes
%
% ========================================================================
% Outputs
% ========================================================================
% M = stiffness matrix

function [M] = fem_assemble(mesh,elem,GQ,shape)

A1 = sparse(mesh.nvN,mesh.nvN);
B1 = sparse(mesh.nvN,mesh.npN);
B2 = B1;
B3 = B1;

for k=1:mesh.nE
    [elem] = fem_calcElem(k,elem,GQ,shape);
    
    A1(mesh.vConnect(k,:), mesh.vConnect(k,:)) = ...
        A1(mesh.vConnect(k,:), mesh.vConnect(k,:)) + elem(k).A1e;
    B1(mesh.vConnect(k,:), mesh.pConnect(k,:)) = ...
        B1(mesh.vConnect(k,:), mesh.pConnect(k,:)) + elem(k).B1e;
    B2(mesh.vConnect(k,:), mesh.pConnect(k,:)) = ...
        B2(mesh.vConnect(k,:), mesh.pConnect(k,:)) + elem(k).B2e;
    B3(mesh.vConnect(k,:), mesh.pConnect(k,:)) = ...
        B3(mesh.vConnect(k,:), mesh.pConnect(k,:)) + elem(k).B3e;
end

I = eye(3);  
A = kron(I,A1);

BT = [B1;B2;B3]; 
B = [B1' B2' B3'];

M = [A BT; B zeros(mesh.npN,mesh.npN)]; % stiffness matrix
end