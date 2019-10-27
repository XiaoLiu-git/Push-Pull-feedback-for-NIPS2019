% Push-Pull demo for real image task
%
%% Original coding by Xiao Liu, Peking University
% $Revision: 1.1.1 $  $Date: 2019-10-27 $
%
% Part of the Push-pull feedback Model version 1.1 for MATLAB.
% See
%"Push-pull Feedback Implements Hierarchical Information Retrieval Efficiently" 
% (NeurIPS 2019), Sec.5 for more details.
% 
%
[raw_instance,Pattern_par] = LoadData(1);%load raw_instance(0,1),Pattern_par(0,1)
[num_gra,num_par,num_chi]=deal(1,2,9);
[N,num_pat]=size(Pattern_par);
Pattern_gra=Gen_highlayer(Pattern_par,num_chi);
wt_par = Hebb_weight( Pattern_par );
wt_gra = Hebb_weight( Pattern_gra );
b1=0.4;
b2=0.17;
tau=5;
B=1;
feedforward= Feedforward(Pattern_gra,Pattern_par,num_gra,num_par);
pfeedback =feedforward';%(xi_b),(xi_c+1)/2
nfeedback=-b1;
Fun=0;
%% inputs of each layers with time
Times=20;
dt=0.01;
time_in=7;
T=0.01; %for sampling

for Num_chi=652;
    m=zeros(9,2,200);
    [I_ext,I_ext2,I_pfb,I_nfb]=deal(zeros(1,Times/dt));
    I_ext(1:Times/dt)=1;
    I_pfb((time_in-2)/dt:(time_in+5)/dt)=1;
    I_nfb((time_in+5)/dt:(time_in+12)/dt)=1;
      
    %% each trail
    [x0]=deal(zeros(N,1));
    [m1,m2,m11]=deal(zeros(1,Times/dt));
    [h1,h2,x1,x2]=deal(zeros(N,1));
    [h11,x11]=deal(zeros(N,1));
    [h1_push,h1_pull,x1_push,x1_pull]=deal(zeros(N,1));
    [n1,n11]=deal(zeros(1,Times/dt));
    [h2_push,h2_pull,x2_push,x2_pull]=deal(zeros(N,1));
    Tn=ceil(T/dt);
    x0=raw_instance(:,Num_chi);  
    for ti=1:Times/T
        for i=1:Tn
            % push-pull feedback
            dh1=1/tau*(1*wt_par*(x1)-h1/B+...
                6*x0*I_ext((ti-1)*Tn+i)+...
                2*pfeedback*(x2)*I_pfb((ti-1)*Tn+i)+...
                1.5*nfeedback*(x2)*I_nfb((ti-1)*Tn+i))*dt;%+...1*(rand(N,1)-0.5)
            h1=h1+dh1;
            x1=0.5*((2/pi)*atan(8*pi*h1)+1);
            % Layer 2
            dh2=1/tau*(2*wt_gra*(x2)-h2/B+...
                1*x0*I_ext((ti-1)*Tn+i)+...
                1*feedforward*(x0))*dt;
            h2=h2+dh2;
            x2=0.5*((2/pi)*atan(8*pi*h2)+1);
            
            
            % Without feedback
            dh11=1/tau*(1*wt_par*(x11)-h11/B+...
                6*x0*I_ext((ti-1)*Tn+i))*dt;%+...           1*(rand(N,1)-0.5)
            h11=h11+dh11;
            x11=0.5*((2/pi)*atan(8*pi*h11)+1);
            % Overlap of Neural activity
            m11((ti-1)*Tn+i)=noverlap(x11,Pattern_par(:,ceil(Num_chi/100)),Fun);
            m1((ti-1)*Tn+i)=noverlap(x1,Pattern_par(:,ceil(Num_chi/100)),Fun);
            m2((ti-1)*Tn+i)=noverlap(x2,Pattern_gra(:,ceil(Num_chi/num_chi/100)),Fun);
            n11((ti-1)*Tn+i)=sum(x11);
            n1((ti-1)*Tn+i)=sum(x1);
  
        end
m(:,:,ti) = m_overlapTi(x1,Pattern_par,num_chi,0);    
    end
    %% save overlap
figure;
Overlaps=plot((dt:dt:Times)./tau,m11,'b-',...
    (dt:dt:Times)./tau,m1,'r-',...
    (dt: dt:Times)./tau,m2,'-');
set(Overlaps,'linewidth',3)
xlabel('Time(\tau)')
ylabel('Retrieval Accuracy')
legend('m_{1,1,1} without fb','m_{1,1,1} push-pull',...
       ' m_{1,1} ')
figure_FontSize=20;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',10),'FontSize',figure_FontSize);
m_overlapTi(x1,Pattern_par,9,1);
end
% figure;
% Overlaps=plot((dt:dt:Times)./tau,n11/4096,'b-',...
%     (dt:dt:Times)./tau,n1/4096,'r-');
% set(Overlaps,'linewidth',3)
% xlabel('Time(\tau)')
% ylabel('Retrieval Accuracy')
% legend('<x> without fb','<x> push-pull')
% figure_FontSize=20;
% set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
% set(findobj('FontSize',10),'FontSize',figure_FontSize);
% figure;
% Fig=plot((-99*dt:dt:Times)./tau,[zeros(1,100),m1],...
%          (-99*dt:dt:Times)./tau,[zeros(1,100),reshape(mean(m(:,1,:)),[1,2000])],...
%          (-99*dt:dt:Times)./tau,[zeros(1,100),reshape(mean(m(:,2,:)),[1,2000])]);
% set(Overlaps,'linewidth',3)
% xlabel('Time(\tau)')
% ylabel('Retrieval Accuracy')
% legend('m^{cat A}','<m>^{cat}','<m>^{dog}')
% figure_FontSize=20;
% set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
% set(findobj('FontSize',10),'FontSize',figure_FontSize);