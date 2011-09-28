function [rval,A]=check_eighistogram(smatfilename,nbins)
A = readSMAT(smatfilename);


if ~exist('nbins','var') || isempty(nbins), nbins=201; end;
[~,bins] = hist([0,2],nbins);

if mod(nbins,2)==0, warning('should call eighistogram with an odd number of bins'); end

v = graph_eigs(A,'normalized');

truebins = hist(v,bins);
eighist = eighistogram(A,nbins);

if any(truebins ~= eighist)
    fprintf('Difference in eigenvalue histogram for %s ("normalized")\n', smatfilename);
    fprintf('  ||true-eigh||_1  = %i\n', norm(truebins - eighist,1));
    fprintf('  ||true-eigh||_oo = %i\n', norm(truebins - eighist,inf));
    rval=false;
else
    fprintf('Passed check on %s ("normalized")\n', smatfilename);
    rval=true;
end

% check adjacency

v = graph_eigs(A,'adjacency');

[eighist,x] = eighistogram(A,nbins,'adjacency');
truebins = hist(v,x);


if any(truebins ~= eighist)
    fprintf('Difference in eigenvalue histogram for %s ("adjacency")\n', smatfilename);
    fprintf('  ||true-eigh||_1  = %i\n', norm(truebins - eighist,1));
    fprintf('  ||true-eigh||_oo = %i\n', norm(truebins - eighist,inf));
    rval=false;
else
    fprintf('Passed check on %s ("adjacency")\n', smatfilename);
    rval=true;
end


% check laplacian

v = graph_eigs(A,'laplacian');

[eighist,x] = eighistogram(A,nbins,'laplacian');
truebins = hist(v,x);


if any(truebins ~= eighist)
    fprintf('Difference in eigenvalue histogram for %s ("laplacian")\n', smatfilename);
    fprintf('  ||true-eigh||_1  = %i\n', norm(truebins - eighist,1));
    fprintf('  ||true-eigh||_oo = %i\n', norm(truebins - eighist,inf));
    rval=false;
else
    fprintf('Passed check on %s ("laplacian")\n', smatfilename);
    rval=true;
end