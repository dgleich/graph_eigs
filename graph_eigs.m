function [v,r,p,Q] = graph_eigs(A,type,novals)

n = size(A,1);

switch type
    case 'adjacency'
        M = A;
    
    case 'laplacian'
        M = diag(sum(A)) - A;
        
    case 'normalized'
        Dhalf = diag(sparse(1./sqrt(sum(A))));
        M = speye(n) - Dhalf*A*Dhalf;
        M = triu(M);
        M = M + triu(M,1)';
        
    case 'modularity'
        d = sum(A)';
        vol = full(sum(d));
        M = A - (1/vol)*(d*d');
        
    otherwise
        error('unknown type');
end

if exist('novals','var') && novals
    v = hess(M);
    return
end    

if nargout > 1
    [V,D] = eig(full(M));
    v = diag(D);
    r = residuals(M,V,v);
    p = sum(V.^4)';
    Q = V;
else
    v = eig(M);
end    

function r = residuals(M,V,d)

r = zeros(size(V,2),1);
for i=1:size(V,2)
    r(i) = norm(M*V(:,i) - d(i)*V(:,i));
end
