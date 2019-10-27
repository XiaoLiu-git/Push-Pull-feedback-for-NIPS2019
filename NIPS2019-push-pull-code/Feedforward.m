function [ feedforward ] = Feedforward(Xi_b,xi_c)
%    Function to generate the feedforward connection
%
% parameter: 
%    Xi_b is the primitive matrix of parents patterns 
%         before the shape changes (size:N,P_a,P_b,P_c)
%    xi_c is the matrix of childern patterns  (size:N,P_a*P_b*P_c)
[~,~,p_b,p_a]=size(Xi_b);
[N,p1]=size(xi_c);
p_c=p1/p_b/p_a;
Xi_b=repmat(Xi_b,[1,p_c,1,1]);
xi_b=reshape(Xi_b,[N,p_a*p_b*p_c]);
xi_b_mean=mean(mean(xi_b,2));xi_b_mean=xi_b_mean*ones(N,p1);
xi_c_mean=mean(mean(xi_c,2));xi_c_mean=xi_c_mean*ones(N,p1);
% xi_b_mean=0;
% xi_c_mean=0;
feedforward=(xi_b-xi_b_mean)*(xi_c-xi_c_mean)';
feedforward=(feedforward-diag(diag(feedforward)))/N;
feedforward=feedforward/norm(feedforward);
end
