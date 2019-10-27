function [ Weight ] = weight( xi )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
[N,p]=size(xi);
Weight=zeros(N,N);
xi_mean=mean(mean(xi,2));xi_mean=xi_mean*ones(N,p);
% xi_mean=0;
Weight=(xi-xi_mean)*(xi-xi_mean)';
Weight=(Weight-diag(diag(Weight)))/N;
Weight=Weight/norm(Weight);
end