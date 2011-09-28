function [neigs,x] = eighistogram(A,nbins,type)
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

if nargin<3, type='normalized'; end

n = size(A,1);

% checks on A
assert(isequal(A,A'));
%assert(all(diag(A)==0));
assert(all(nonzeros(A)==1));
assert(size(A,1)==size(A,2));

% define what zero should be for finding the histogram bins
% perhaps make this an option at some point.
logzero = log10(100*eps);

switch type
    case 'normalized'
        Dhalf = diag(sparse(1./sqrt(sum(A))));
        M = speye(n) - Dhalf*A*Dhalf;
        M = triu(M);
        M = M + triu(M,1)';
        % compute the histogram bins
        [ignore,x] = hist([0,2],nbins); %#ok<ASGLU>

    case 'adjacency'
        M = A;
        dmax = full(max(sum(A)));
        assert(mod(nbins,2)==1); 
        nhalf = (nbins+1)/2;
        xhalf = logspace(logzero,log10(dmax), nhalf);
        x = [-fliplr(xhalf(2:end)) xhalf];
        
    case 'laplacian'
        dmax = full(max(sum(A)));
        M = diag(sum(A)) - A;
        % use a slightly higher tolerance for zero in this 
        % call.  Empirically, eigenvalues of less than 1e-12 
        % correlate with disconnected components.
        x = [0 logspace(logzero+1.5,log10(2*dmax),nbins-1)];
        
    otherwise
        error('graph_eigs:eighistogram','unknown type "%s"', type);
end

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
% Try two at an algorithm
[nzi,nzj,nzv] = find(D);
pos = 0;
neg = 0;
zero = 0;
d = zeros(size(D,1),1); % flag for a block-diagonal we saw
for nz=1:length(nzi)
    % get the element
    i = nzi(nz);
    j = nzj(nz);
    v = nzv(nz);
    if i~=j,
        d(i)=true; 
        d(j)=true;
    end
end
for nz=1:length(nzi)
    % get the element
    i = nzi(nz);
    j = nzj(nz);
    v = nzv(nz);
    if ~d(i) && ~d(j)
        % this row/col pair isn't in any off-diagonal block
        assert(i==j);
        if v < eps(1)
            neg = neg + 1;
        elseif v > eps(1)
            pos = pos + 1;
        else
            zero = zero + 1;
        end
    end
end
% add in additional zero eigenvalues
zero = zero + (sum(~d)-pos-neg-zero);

B = zeros(2,2);
bfull = false;
bstart = 0;
bend = 0;
for nz=1:length(nzi)
    i = nzi(nz);
    j = nzj(nz);
    v = nzv(nz);
    
    if d(i) || d(j)
        assert(d(i) && d(j)); % both of them have to be flagged
        if i>bend || j>bend
            % we've finished the block, compute it's eigenvalues
            if bfull
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
            end
            
            bstart=min(i,j);
            bend=bstart+1;
            B(i-bstart+1,j-bstart+1) = v;
            
            bfull=true;
        else
            % we're inputting stuff from the block
            B(i-bstart+1,j-bstart+1) = v;
        end
    end
end

% final update
if bfull
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
end
%[pos, neg, zero, pos+neg+zero, size(D,1)]
assert(pos+neg+zero == size(D,1));    

