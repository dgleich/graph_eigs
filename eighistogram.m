function [neigs,x] = eighistogram(A,nbins)
% EIGHISTOGRAM Compute eigenvalue density using the ldl factorization
%
% This function only works on the normalized Laplacian matrix at the
% moment.  
% 
% It uses the inertial preserving properties of the LDL' factorization to
% count the number of eigenvalues in a region given a set of shifts.
%
% [n,x] = eighistogram(A,nbins) will form the normalized Laplacian matrix
% for the graph with adjacency matrix A, then it will compute nbins-1 LDL
% factorizations to count the eigenvalues.
% 
% Example:
%   A = readSMAT('test/Caltech36.smat');
%   [n,x] = eighistogram(A,35);
%   v = graph_eigs(A,'normalized');
%   [nv,xv] = hist(v,35);
%   plot(x,n,'o-',xv,nv,'.--'); % yes, they are right on top of each other

n = size(A,1);

Dhalf = diag(sparse(1./sqrt(sum(A))));
M = speye(n) - Dhalf*A*Dhalf;
M = triu(M);
M = M + triu(M,1)';

% compute the histogram
[ignore,x] = hist([0,2],nbins);
binedges = x(1:end-1) + diff(x)/2;

neigs = zeros(1,length(x));

cureigs = 0;

if n > 30000
    verbose = true;
else 
    verbose = false;
end

for i=1:length(binedges)
    
    be = binedges(i);
    T = M - be*speye(n);
    if verbose
        fprintf('%4i/%4i computing ldl for shift %6f ... ', i,length(binedges), be);
    end
    t0 = tic;
    [L,D,p,s] = ldl(T);
    dt = toc(t0);
    if verbose
        fprintf('%.1f secs.\n', dt);
    end
    [pos,neg,zero] = inertia(D);
    neigs(i) = neg + zero - cureigs;
    cureigs = neigs(i) + cureigs;
end
neigs(end) = n - cureigs; % all that remain are at the end

function [pos,neg,zero] = inertia(D)
% INERTIA Compute the inertia of a diagonal or 2-by-2 block diagonal pair

% Todo: put in normalized zero check
[nzi,nzj,nzv] = find(D);
pos = 0;
neg = 0;
zero = 0;
d = zeros(size(D,1),1); % flag for a diagonal or block-diagonal we saw

% the idea in the next part is to look over all the non-zeros to compute
% the inertia.  we'll do this by following the algorithm:
%   for each non-zero
%     if in the same column or off-diagonal, 
%       then start 2x2 block
%     else
%        it was 

% initialize our running structure

% twobytwo if we are in a full 2x2 block
% twobytwooff if we are in a 2x2 off-diagonal block

% todo handle empty d

twobytwo = false;
twobytwoelem = 1;
twobytwooff = nzi(1)~=nzj(1);
d(nzj(1)) = 1;
lastelemnew = true; % true if the last element needs to be handled still
lastelem = nzv(1);
lastcol = nzj(1); % the last column we saw.
B = zeros(2); % a diagonal block, as such.

for nz=2:length(nzi)
    % get the element
    i = nzi(nz);
    j = nzj(nz);
    v = nzv(nz);
    d(j) = 1;
    if twobytwooff
        assert(i ~= j);
        if abs(lastelem) < eps(1)
            zero = zero + 2;
        else
            pos = pos + 1;
            neg = neg + 1;
        end
        twobytwooff = false;
        lastelemnew = false;
    elseif twobytwo
        assert(twobytwoelem > 2); % only detected after the second element
        B(twobytwoelem) = v;
        twobytwoelem = twobytwoelem + 1;
        if twobytwoelem > 4
            % we are done!
            eigs = eig(B);
            B(:) = 0; % reset B
            
            % update counts
            if eigs(1) < eps(1)
                neg = neg + 1;
            elseif eigs(1) > eps(1)
                pos = pos + 1;
            else
                zero = zero + 1;
            end
            
            if eigs(2) < eps(1)
                neg = neg + 1;
            elseif eigs(2) > eps(1)
                pos = pos + 1;
            else
                zero = zero + 1;
            end
            
            twobytwo = false;
            twobytwoelem = 1;
            lastelemnew = false;
        end
    else
        % cases:
        % 1. we see a new column
        if j > lastcol
            % if the last element hasn't already been processed, hit it now
            if lastelemnew
                if lastelem < eps(1)
                    neg = neg + 1;
                elseif lastelem > eps(1)
                    pos = pos + 1;
                else
                    zero = zero + 1;
                end
            end
            
            % otherwise, if we find a diagonal element
            if i==j
                % then we don't know what to do, it oculd be the start
                % of a block
                lastelemnew = true;
            else
                % ahah!  A new column _and_ an off-diagonal, we have
                % an off-diagonal 2x2 block
                twobytwooff = true;
                lastelemnew = false;
            end
        else
            % okay, here we are starting a 2x2 block
            twobytwo = true;
            B(1) = lastelem;
            B(2) = v;
            twobytwoelem = 3;
            lastelemnew = false;
        end
    end
    lastelem = v;
    lastcol = j;
end
% handle the final element
if lastelemnew
    if lastelem < eps(1)
        neg = neg + 1;
    elseif lastelem > eps(1)
        pos = pos + 1;
    else
        zero = zero + 1;
    end
end
zero = zero + (size(D,1) - sum(d));
%[pos, neg, zero, size(D,1)]
assert(pos+neg+zero == size(D,1));    