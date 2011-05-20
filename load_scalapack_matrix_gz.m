function A = load_scalapack_matrix_gz(filename)
% Load a matrix written by pdlawrite, then compressed with gzip
% A = load_scalpack_matrix(filename)
% the file format is 
% <m> <n>
% <val>@(<m>*<n>)

% Create Java input stream from the gzipped filename.
fileInStream = [];
try
   fileInStream = java.io.FileInputStream(java.io.File(filename));
catch exception
   % Unable to access the gzipped file.
   if ~isempty(fileInStream)
     fileInStream.close;
   end
   eid = sprintf('MATLAB:%s:javaOpenError',mfilename);
   error(eid,'Could not open file "%s" for reading.',filename);
end

% Create a Java GZPIP input stream from the file input stream.
try
   gzipInStream = java.util.zip.GZIPInputStream( fileInStream );
catch exception
   % The file is not in gzip format.
   if ~isempty(fileInStream)
     fileInStream.close;
   end
   eid = sprintf('MATLAB:%s:notGzipFormat',mfilename);
   error(eid,'File "%s" is not in GZIP format.',filename);
end

gzipChannel = java.nio.channels.Channels.newChannel(gzipInStream);
bufsize = 16*1024*1024; % 16MB buffer
B = java.nio.ByteBuffer.allocate(bufsize);

% read a chunk
B.clear();
len = gzipChannel.read(B);
str = char(B.array()');
filedone = false;
if len < length(str)
    str = str(1:len);
    filedone = true;
end
% find the 
if ~filedone
    lastnl = find_last_newline(str);
    strextra = str(lastnl+1:end);
    str = str(1:lastnl);
else 
    strextra = [];
end

% read the header
[s,pos] = textscan(str, '%d',2);
str = str(pos+1:end); % truncate the string to what was used.

m = s{1}(1);
n = s{1}(2);

A = zeros(m, n);
curind = 1;

while 1
    
    [buf,pos] = textscan(str,'%f');
    buf = buf{1};
    A(curind:curind+length(buf)-1) = buf;
    curind = curind + length(buf);
    
    str = str(pos+1:end);
    assert(isempty(str)); % we should have used all of the string
    
    if filedone
        % if filedone was set, we are done!
        break
    end
    
    % read the next file chunk
    B.clear();
    len = gzipChannel.read(B);
    str = char(B.array()');
    filedone = false;
    if len < length(str)
        str = str(1:len);
        filedone = true;
    end
    % add in the old extra
    if ~isempty(strextra)
        str = [strextra str];
    end
    if ~filedone
        % find the last newline
        lastnl = find_last_newline(str);
        strextra = str(lastnl+1:end);
        str = str(1:lastnl);
    else 
        strextra = [];
    end
end

if ~isempty(fileInStream)
 fileInStream.close;
end


function ind=find_last_newline(str)

ind = [];
for i=length(str):-1:1
    if str(i) == 10 || str(i) == 13 % these are the newline characters
        ind = i;
        break
    end
end

% s = textscan(fid, '%d',2);
% m = s{1}(1);
% n = s{1}(2);
% A = textscan(fid, '%f', m*n,'CollectOutput',1);
% A = reshape(A{1}, m, n);
% fclose(fid);