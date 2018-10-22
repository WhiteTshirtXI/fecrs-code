function dY = rates_fem_rsm(t,Y,Minv,mesh,bnd,epsilon,save_data,nSave)
% attempt #1 of encoding FEM-RSM method such that an adaptive ode solver
% such as ode15s can be used (rather than fixed-dt forward Euler method).

% This code should take all input required to compute the FEM solution, producing a
% vector of velocities at bead location.
% This velocity data will then be fed through ode15s to compute bead locations.

% Input vector Y has the Q+1 bead locations.

beads           = Y;                    % column vector formed of bead locations (x1,x2,x3)'
Q               = (length(Y)-3)/3;      % number of segments
N               = Q+1;                  % number of beads/joints
beads_comp(:,1) = beads(1:N);           % component-wise form of bead locations
beads_comp(:,2) = beads(N+1:2*N);
beads_comp(:,3) = beads(2*N+1:end);

% solve blm for bead and boundary velocities
[Vf,Vb] = rates_rsm_fecrs(beads(:),epsilon,bnd);

% calculate load vector for fem
[b] = fem_calcLoadVec(mesh,Vb);

us = Vf(1:(Q+1));
vs = Vf(Q+2:(2*Q+2));
ws = Vf(2*Q+3:(3*Q+3));

% solve system
c = Minv*b;

% velocity solutions at nodes of fem mesh
U(:) = c(1:mesh.nvN,1);                 % u solution
V(:) = c(1+mesh.nvN:2*mesh.nvN,1);      % v solution
W(:) = c(2*mesh.nvN+1:3*mesh.nvN,1);    % w solution
%P(:) = c((3*mesh.nvN+1):end);          % p solution

% evaluate fem solution at all beads
[U_bead,V_bead,W_bead] = fem_evalVel(mesh,beads_comp(:,1),beads_comp(:,2),beads_comp(:,3)...
                                    ,U(:),V(:),W(:));

% evaluate final velocity at beads (u = U + us)
for i=1:N
    u_bead(i) = U_bead(i) + us(i);
    v_bead(i) = V_bead(i) + vs(i);
    w_bead(i) = W_bead(i) + ws(i);
end

dY = [u_bead(:);v_bead(:);w_bead(:)]; 

% s(t+1)  = calc_arclength(X(:,1:2,t)'); % calculate arclength
% perccount(t,Nt); % displays percentage complete

% if (save_data == 1) 
%     if (mod(t,nSave) == 0) % save every nSave timesteps
%         Xsave(:,:,count) = X(:,:,t);
%         save(filename,'Xsave','tvec','nBeads','u_bead','v_bead','w_bead','s'); 
%         count = count + 1;
%     end
% end

end % function