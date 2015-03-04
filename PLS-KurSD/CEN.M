function y=cen(x,m) 
if nargin < 2, m=mean(x); end 
u=ones(size(x)); 
y=x-u(:,1)*m; 
