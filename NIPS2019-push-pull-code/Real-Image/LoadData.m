function [raw_instance,Pattern_par] = LoadData(x)
% Function to loading data
% parameter: 
%     ifx==1 then laod data, else not
if x==1
load('simpleCD.mat')% neural activities filtered by VGG
% normalization 
raw_instance=log(1+double(fc2s'));
Mean_fc2s=ones(4096,1)*mean(raw_instance,1);
Std_fc2s=ones(4096,1)*std(raw_instance,0,1);
Std_fc2s(Std_fc2s==0)=1;
raw_instance=(raw_instance-Mean_fc2s)./(2*Std_fc2s)+0.5;%
raw_instance(raw_instance<0)=0;
load('Pattern_par.mat')%Pattern_par
Mean_fc2s=ones(4096,1)*mean(raw_instance,1);
raw_instance=sign(raw_instance-Mean_fc2s)*0.5+0.5;
Pattern_par=sign(Pattern_par-ones(4096,1)*mean(Pattern_par))/2+0.5;
end
end

