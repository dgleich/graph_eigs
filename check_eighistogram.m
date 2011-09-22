function [rval,truebins,eighist]=check_eighistogram(smatfilename,nbins)
A = readSMAT(smatfilename);


if ~exist('nbins','var') || isempty(nbins), nbins=201; end;
[~,bins] = hist([0,2],nbins);

if mod(nbins,2)==0, warning('should call eighistogram with an odd number of bins'); end

v = graph_eigs(A,'normalized');

truebins = hist(v,bins);
eighist = eighistogram(A,nbins);

if any(truebins ~= eighist)
    fprintf('Difference in eigenvalue histogram for %s\n', smatfilename);
    fprintf('  ||true-eigh||_1  = %i\n', norm(truebins - eighist,1));
    fprintf('  ||true-eigh||_oo = %i\n', norm(truebins - eighist,inf));
    rval=false;
else
    fprintf('Passed check on %s\n', smatfilename);
    rval=true;
end