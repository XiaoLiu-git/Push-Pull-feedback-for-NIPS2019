function [ pattern_par ] = Gen_highlayer( pattern_chi,num_chi )
%  Generate patter in higer layer
%
% parameter: 
%     pattern_chi: pattern in the lower layer
%     num_chi: the number of childern for each parent
[N,num_pat]=size(pattern_chi);
pattern_par=reshape(pattern_chi,[N,num_chi,num_pat/num_chi ]);
pattern_par=reshape(mean(pattern_par,2),[N,num_pat/num_chi]);
pattern_par=sign(pattern_par-mean(mean(pattern_par)))/2+0.5;
%% Change all the elements from 0.5 to 0
I= pattern_par==0.5; 
pattern_par(I)=0;
end

