% Generates nodes and weights for 2x2x2 Gaussian quadrature rule
%
% [GQ] = fem_setupQuadrature()
%
% =======================================================================
% Outputs
% =======================================================================
% GQ.node = nodes for Gaussian quadrature rule
% GQ.weight = weights for Gaussian quadrature rule

function [GQ] = fem_setupQuadrature()

% nodes
GQ.node(1,1) = -1/sqrt(3);
GQ.node(2,1) = 1/sqrt(3);
GQ.node(3,1) = 1/sqrt(3);
GQ.node(4,1) = -1/sqrt(3);
GQ.node(5,1) = -1/sqrt(3);
GQ.node(6,1) = 1/sqrt(3);
GQ.node(7,1) = 1/sqrt(3);
GQ.node(8,1) = -1/sqrt(3);

GQ.node(1,2) = -1/sqrt(3);
GQ.node(2,2) = -1/sqrt(3);
GQ.node(3,2) = 1/sqrt(3);
GQ.node(4,2) = 1/sqrt(3);
GQ.node(5,2) = -1/sqrt(3);
GQ.node(6,2) = -1/sqrt(3);
GQ.node(7,2) = 1/sqrt(3);
GQ.node(8,2) = 1/sqrt(3);

GQ.node(1,3) = -1/sqrt(3);
GQ.node(2,3) = -1/sqrt(3);
GQ.node(3,3) = -1/sqrt(3);
GQ.node(4,3) = -1/sqrt(3);
GQ.node(5,3) = 1/sqrt(3);
GQ.node(6,3) = 1/sqrt(3);
GQ.node(7,3) = 1/sqrt(3);
GQ.node(8,3) = 1/sqrt(3);

% weights
GQ.weight(1) = 1;
GQ.weight(2) = 1;
GQ.weight(3) = 1;
GQ.weight(4) = 1;
GQ.weight(5) = 1;
GQ.weight(6) = 1;
GQ.weight(7) = 1;
GQ.weight(8) = 1;

end