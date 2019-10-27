function [xi_a,xi_b,xi_c] = HirPatterns(N,p_a,p_b,p_c,b1,b2)
%   HirPatterns will return a set of hierarchical patterns with p_b classes
%   and p_c patterns in each class.
%   each pattern has 1*N sites. 
%  
% parameter: 
%            each pattern has 1*N sites 
%            b1=0~1 is the correlation bewteen two patterns in the same class.
%            b2=0~1 is the correlation of two class patterns

%% Generate {-1,1} patterns %%
%%xi_a%%
xi_a=binornd(1,0.5,[N,1,1,p_a])*2-1;

%% xi_b %%
xi_b=rand(N,1,p_b,p_a);
xi_b=xi_b-(1-b2*repmat(xi_a,[1 1 p_b 1]))/2>0;
xi_b=double(xi_b)*2-1;
xi_a=sign(sum(xi_b,3));

%% xi_c %%
xi_c=rand(N,p_c,p_b,p_a);
xi_c=xi_c-(1-b1*repmat(xi_b,[1 p_c 1 1]))/2>0;
xi_c=double(xi_c)*2-1;
xi_b=sign(sum(xi_c,2));
xi_a=sign(sum(xi_b,3));
%% Change all the elements from 0 to 1
I=find(xi_a==0);
xi_a(I)=-1;
I=find(xi_b==0);
xi_b(I)=-1;
I=find(xi_c==0);
xi_c(I)=-1;
%% Align patterns 
xi_a=(xi_a+1)/2;
xi_b=(xi_b+1)/2;
xi_c=(xi_c+1)/2;
end
