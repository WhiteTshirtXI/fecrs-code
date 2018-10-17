% Calculates Jacobian matrix and its determinant for each element. Also
% calculates the derivatives of the velocity and pressure shape functions
% w.r.t x,y and z at GQ points for each element

% [elem] = fem_calcJacob(mesh,dshape)
%
% ========================================================================
% Inputs
% ========================================================================
% mesh = structure containing mesh parameters
% dshape = structure containing derivatives of both velocity and pressure
% shape functions evaluated at GQ nodes
%
% ========================================================================
% Outputs
% ========================================================================
% elem.detJacob.p = determinant of pressure Jacobain for each element
% elem.detJacob.v = determinant of velocity Jacobain for each element
% elem.gDS.p = derivatives of pressure shape functions w.r.t x and y at GQ
% nodes
% elem.gDS.v = derivatives of velocity shape functions w.r.t x and y at GQ
% nodes

function [elem] = fem_calcJacob(mesh,dshape)

for k=1:mesh.nE
    
    xp = mesh.pCoords(mesh.pConnect(k,:),1);
    yp = mesh.pCoords(mesh.pConnect(k,:),2);
    zp = mesh.pCoords(mesh.pConnect(k,:),3);
    
    xv = mesh.vCoords(mesh.vConnect(k,:),1);
    yv = mesh.vCoords(mesh.vConnect(k,:),2);
    zv = mesh.vCoords(mesh.vConnect(k,:),3);
    
    e_coords.p = [xp yp zp];
    e_coords.v= [xv yv zv];
    
    for j=1:8
        % jacobians
        Jacob.p(:,:) = dshape.p(:,:,j)*e_coords.p(:,:);
        Jacob.v(:,:) = dshape.v(:,:,j)*e_coords.v(:,:);
        
        % determinants of Jacobians
        elem(k).detJacob.p(j) = det(Jacob.p);
        elem(k).detJacob.v(j) = det(Jacob.v);
        
        % shape function derivatives w.r.t x and y
        elem(k).gDS.p(:,:,j) = inv(Jacob.p(:,:))*dshape.p(:,:,j);
        elem(k).gDS.v(:,:,j) = inv(Jacob.v(:,:))*dshape.v(:,:,j);
    end
end

end