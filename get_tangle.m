function theta = get_tangle(X)
% calculates the tangent angle of a given vector X, or horizontal
% concatenation of vectors X.
% input vectors should be 3D ie X = [x1; x2; x3].

for n = 1:size(X,2)
    nX            = norm(X(:,n));
    theta(n,1)    = acos( dot(X(:,n),[1;0;0])/nX );
end

end % function