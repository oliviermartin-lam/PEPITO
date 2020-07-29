function [result, rr,cum] = radial(matrix, m, n)
% Find radial average of input matrix
% Input matrix must be square.

if nargin < 2
    [m,n] = size(matrix);
    centerm = ceil(m/2+1); %matrix is square , so use m or n
    centern = centerm;
else
	centerm = n;
    centern = m;
    [m,n] = size(matrix);
end

[i,j] = meshgrid(1:n,1:m);
i = i(:); j = j(:);

dist = (sqrt((i-centerm).^2 + (j-centern).^2));
[dist,y] = sort(dist);

% this is because I need to find how many matrix values are located at
%each distance
% so I decided to use 'hist' function to find that
hh = hist(dist,max(dist)+1);

vec = matrix(:);
vec = vec(y); % sort the same way as dist

ini = 1;

result(1:floor(round(max(dist))))=0;
cum(1:floor(round(max(dist))))=0;
for k = 1:floor(max(dist)) %maybe this loop still can be optmized

   index = ini:ini+hh(k+1)-1; % this indexing is working perfectly
   result(k) = mean(vec(index));
   cum(k) = sum(vec(index));
   ini = max(index)+1;

end

result = result(1:floor(min(size(matrix))/2)-1);
result = [matrix(floor(centern),floor(centerm)) result];
cum = cum(1:floor(min(size(matrix))/2)-1);
cum = [matrix(floor(centern),floor(centerm)) cum];

rr = unique(dist);
rr = rr(1:length(result));