function [ m ] = moverlap(x,xi,method)
%UNTITLED3 Summary of this function goes here
%   method defines the function for calculating the overlap,when
%   method=0, we use sign(x), else we use sign atan(method*x)
N=length(x);
if method==0
m=1/N*(sign(x*2-1)'*sign(xi*2-1));
else
m=1/N*((2/pi)*atan(method*(x*2-1))'*sign(xi*2-1));

end
