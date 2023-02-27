% p = ncc(x,y)
%
% Calculate the normalized cross correlation between two matrices x and y
function p = ncc(x,y)

assert(all(size(x)==size(y)),'Size of x and y must be the same')

x=x(:);x=x-mean(x);
y=y(:);y=y-mean(y);
p=x'*y/sqrt((x'*x)*(y'*y));
