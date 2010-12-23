function A = load_scalapack_matrix(filename)
% Load a matrix written by pdlawrite
% A = load_scalpack_matrix(filename)
% the file format is 
% <m> <n>
% <val>@(<m>*<n>)

fid = fopen(filename);
if fid == -1
    error('Cannot open file');
end
s = textscan(fid, '%d',2);
m = s{1}(1);
n = s{1}(2);
A = textscan(fid, '%f', m*n,'CollectOutput',1);
A = reshape(A{1}, m, n);
fclose(fid);