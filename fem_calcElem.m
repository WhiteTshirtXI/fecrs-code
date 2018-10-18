% Calculates element stiffness matrices to be used in the construction
% of the whole stiffness matrix
%
% [elem] = fem_calcElem(e,elem,GQ,shape)
%
% ========================================================================
% Inputs
% ========================================================================
% e = element number
% elem = structure consisting of jacobian matrices e.t.c required to map
% from reference element
% GQ = structure containing required Gaussian quadrature nodes and weights
% shape = structure containing velocity and pressure shape functions
% evaluated at GQ nodes
%
% ========================================================================
% Outputs
% ========================================================================
% elem = structure containing element stiffness matrices

function [elem] = fem_calcElem(e,elem,GQ,shape)

elem(e).A1e = zeros(27,27);
elem(e).B1e = zeros(8,27);
elem(e).B2e = zeros(8,27);
elem(e).B3e = zeros(8,27);

for k=1:8
    elem(e).A1e(:,:) = elem(e).A1e(:,:) + ...
        (elem(e).gDS.v(1,:,k) .* elem(e).gDS.v(1,:,k)' + ...
        elem(e).gDS.v(2,:,k) .* elem(e).gDS.v(2,:,k)' + ...
        elem(e).gDS.v(3,:,k) .* elem(e).gDS.v(3,:,k)') ...
        .* elem(e).detJacob.v(k) .* GQ.weight(k);
    
    elem(e).B1e(:,:) = elem(e).B1e(:,:) + ...
        (elem(e).gDS.v(1,:,k).*shape.p(:,k)) ...
        .* elem(e).detJacob.v(k) .* GQ.weight(k);
    
    elem(e).B2e(:,:) = elem(e).B2e(:,:) + ...
        (elem(e).gDS.v(2,:,k).*shape.p(:,k)) ...
        .* elem(e).detJacob.p(k) .* GQ.weight(k);
    
    elem(e).B3e(:,:) = elem(e).B3e(:,:) + ...
        (elem(e).gDS.v(3,:,k).*shape.p(:,k)) ...
        .* elem(e).detJacob.p(k) .* GQ.weight(k);
    
end

elem(e).B1e = elem(e).B1e(:,:)';
elem(e).B2e = elem(e).B2e(:,:)';
elem(e).B3e = elem(e).B3e(:,:)';
end