% Push-Pull demo for continuous Hopfield model
%
%% Original coding by Xiao Liu, Peking University
% $Revision: 1.1.1 $  $Date: 2019-10-27 $
%
% Part of the Push-pull feedback Model version 1.1 for MATLAB.
% See
%"Push-pull Feedback Implements Hierarchical Information Retrieval Efficiently" 
% (NeurIPS 2019), Sec.4 for more details.
% 
%
clear all
%% Parameter setting
a_ext=1;  % The strength of the external input
N=2000;   % The number of neurons in a layer
p_a=2;    % The number of grandparents pattern --> P_alpha
p_b=4;    % The number of parents pattern --> P_beta
p_c=25;   % The number of childern pattern --> P_gamma
b1=0.15;  % The strength of correlation between childern patterns, 0<b1<1.
b2=0.1;   % The strength of correlation between parents patterns, 0<b2<1.
tau=5;    % The time constant.
B=1;      % The strength of the recovre rate

for trail = 1:10;
    [xi_a,Xi_b,xi_c] = HirPatterns(N,p_a,p_b,p_c,b1,b2);
    %Xi_b is the primitive matrix before the shape changes
    xi_a=reshape(xi_a,[N,p_a]);
    xi_b=reshape(Xi_b,[N,p_a*p_b]);
    xi_c=reshape(xi_c,[N,p_a*p_b*p_c]);
    
    %% interaction between layers
    W2= weight(xi_b);
    W1= weight(xi_c);
    % feedforward= c_ff*feedforward((xi_b+1)/2,(xi_c+1)/2);
    % pfeedback = c_fb*p_feedback((xi_b+1)/2,(xi_c+1)/2);
    feedforward= Feedforward(Xi_b,xi_c);
    pfeedback =feedforward';%(xi_b),(xi_c+1)/2
    nfeedback=-b1;
    % feedforward=b1^3*(p_c-1);
    % pfeedback = b1^3*(p_c-1);
    clear Xi_b
    
    %% inputs of each layers with time
    time_in=7;
    Time_total=35;dt=0.01;
    X_00=33;
    [x0]=deal(zeros(N,1));
    % x0(1:1*N,1)=xi_c(1:1*N,X_00);
    % x0(0.9*N:1*N,1)=xi_c(0.9*N:1*N,X_00+2*p_c);
    [I_ext,I_ext2,I_pfb,I_nfb]=deal(zeros(1,Time_total/dt));
    I_ext(5/dt:Time_total/dt)=1;
    I_pfb((time_in-2+5)/dt:(time_in+3+5)/dt)=1;
    I_nfb((time_in+3+5)/dt:(time_in+8+5)/dt)=1;
    % n_x2=randperm(length(x1),0.3*N);x2(n_x2,1)=-x1(n_x2,1);
    % x2=x0;
    lambda1=0.2;
    lambda2=0.1;
    
    %% each trail
    [m1,m2,m11,X1_mean,X11_mean]=deal(zeros(1,Time_total/dt));
    x0=xi_c(:,X_00);
    n1_x0=randperm(length(x0),fix(lambda1*N));
    x0(n1_x0,1)=xi_c(n1_x0,X_00+2);
    n_x0=randperm(length(x0),fix(lambda2*N));
    x0(n_x0,1)=xi_c(n_x0,X_00+1*p_c);
    
    [h1,h2,x1,x2]=deal(zeros(N,1));
    [h11,x11]=deal(zeros(N,1));
    for i=5/dt:Time_total/dt
        dh1=1/tau*(1*W1*(x1)-h1/B+...
            a_ext*x0*I_ext(i)+...
            1*pfeedback*(x2)*I_pfb(i)+...
            10*nfeedback*(x2)*I_nfb(i)+...
            0.1*(rand(N,1)-0.5))*dt;
        h1=h1+dh1;
        x1=0.5*((2/pi)*atan(8*pi*h1)+1);
        
        dh11=1/tau*(1*W1*(x11)-h11/B+...
            a_ext*x0*I_ext(i)+...
            0.1*(rand(N,1)-0.5))*dt;
        h11=h11+dh11;
        x11=0.5*((2/pi)*atan(8*pi*h11)+1);
        
        dh2=1/tau*(2*W2*(x2)-h2/B+...
            2*feedforward*(x1)+...
            0*x0*I_ext(i))*dt;%+0.01*(rand(N,1)-0.5)
        h2=h2+dh2;
        x2=0.5*((2/pi)*atan(8*pi*h2)+1);
        
        m11(i)=1/N*(sign(x11*2-1)'*sign(xi_c(:,X_00)*2-1));
        m1(i)=1/N*(sign(x1*2-1)'*sign(xi_c(:,X_00)*2-1));
        m2(i)=1/N*(sign(x2*2-1)'*sign(xi_b(:,ceil(X_00/p_c))*2-1));
        X1_mean(i)=mean(x1,1);
        X11_mean(i)=mean(x11,1);
    end
    figure();
    subplot(1,2,1);
    Overlaps=plot((dt:dt:Time_total)/tau-1,m11,'b',...
        (dt:dt:Time_total)/tau-1,m1,'r',...
        (dt:dt:Time_total)/tau-1,m2,'y');
    set(Overlaps,'linewidth',3)
    xlabel('Time(\tau)')
    ylabel('Retrieval Accuracy')
    legend('m^1 without fb',' m^1 ',' m^2 ')
    figure_FontSize=20;
    set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
    set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
    set(findobj('FontSize',10),'FontSize',figure_FontSize);
%     saveas(gcf,[Path,'/','MisingCurve',num2str(trial),'.eps'],'psc2');
    subplot(1,2,2);
    Overlaps=plot((dt:dt:Time_total)/tau-1,X11_mean,'b',...
        (dt:dt:Time_total)/tau-1,X1_mean,'r');
    set(Overlaps,'linewidth',3)
    ylabel('Activity ⟨x⟩')
    xlabel('Time(\tau)')
    legend('<x^1> without fb',' <x^1> ')
    figure_FontSize=20;
    set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
    set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
    set(findobj('FontSize',10),'FontSize',figure_FontSize);
%     saveas(gcf,[Path,'/','MisingCurve',num2str(trial),'.eps'],'psc2');
end
