function check_commute_times(fn)
% CHECK_COMMUTE_TIMES Check the commute times computed by graph_eigs
% fn - the filename of the smat file, the commute times are stored in 
%   fn.laplacian.ctimes

smatdata = load(fn);
edges = smatdata(2:end,:);
G = spones(sparse(edges(:,1)+1,edges(:,2)+1,1,smatdata(1,1),smatdata(1,2)));
L = diag(sum(G,2))-G;
C = pinv(full(L));
d = diag(C);
n = size(G,1);
for i=1:n
    for j=1:n
        C(i,j) = d(i) + d(j) - 2*C(i,j);
    end
end

C2 = load_scalapack_matrix([fn '.laplacian.ctimes']);

diffs = abs(C(:)-C2(:))./max(C(:),1);
errs = sum(diffs>n*eps(1));
if errs/n^2>0.02
    fprintf('Error, %.1f \% of commute times have high error\n', 100*errs/n^2);
else
    fprintf('Okay!');
end
