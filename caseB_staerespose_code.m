clc
clear
close all
tic

n=5; 
D=1; M=10; Tt=0.3; Tg=0.1; R=0.05; bbeta=21; rhoe=1/21; kappae=1.0; Te=0.1; alpha1=0.5; alpha2=0.5;

Kp1= 0.05;    
Ki1=0.1;
Kp2=0.05; 
Ki2=0.1;


A1=[-D/M,   1/M,     0,   1/M, 0;
    0, -1/Tt,  1/Tt,   0,  0;
    -1/(R*Tg),     0, -1/Tg,  0,  0;
    -(rhoe*kappae)/Te,  0,  0, -1/Te, 0;
    bbeta,     0,     0,   0,   0];
A2=[-D/M,   1/M,     0,   1/M, 0;
    0, -1/Tt,  0.5/Tt,   0,  0;
    -1/(R*Tg),     0, -0.5/Tg,  0,  0;
    -(rhoe*kappae)/Te,  0,  0, -1/Te, 0;
    bbeta,     0,     0,   0,   0];
Ad1=[0,   0,     0,   0,  0;
    0,   0,     0,   0,  0;
    -(alpha1*Kp1*bbeta)/Tg,   0,  0,  0, -(alpha1*Ki1)/Tg;
    -(alpha2*kappae*Kp1*bbeta)/Te,   0,  0,  0, -(alpha2*kappae*Ki1)/Te;
    0,   0,     0,   0,  0];
Ad2=[0,   0,     0,   0,  0;
    0,   0,     0,   0,  0;
    -(alpha1*Kp2*bbeta)/Tg,   0,  0,  0, -(alpha1*Ki2)/Tg;
    -(alpha2*kappae*Kp2*bbeta)/Te,   0,  0,  0, -(alpha2*kappae*Ki2)/Te;
    0,   0,     0,   0,  0];

C=[bbeta,0,0,0,0;
    0,0,0,0,1];


tau1=@(t)0.5*sin(t)*sin(t)+10;

d=0.02; 
t_begin=0;
t_end=80;
ta=-50:d:t_begin; Ta=length(ta);
tb=d:d:t_end;     t=[ta,tb];
T=length(t);

% X=zeros(n,T);
X1=zeros(n,T);
% X2=zeros(n,T);
X0=[0.2; 0.1; -0.3; -0.5; 0.1]; 

for i=1:Ta
    % X(:,i)=X0;
    X1(:,i)=X0;
    %  X2(:,i)=X0;
    % phi=0;
end

for p=Ta:T
    Tau1=fix((tau1((p-Ta)*d))/d);
    X1(:,p+1)=X1(:,p)+ 0.5*(d*A1*X1(:,p)+d*Ad1*X1(:,p-Tau1))+ 0.5*(d*A2*X1(:,p)+d*Ad2*X1(:,p-Tau1));%%系统表达式
   
end

figure;

plot(tb,X1(1,Ta:T-1));hold on;
plot(tb,X1(2,Ta:T-1));hold on;
plot(tb,X1(3,Ta:T-1));hold on;
plot(tb,X1(4,Ta:T-1));hold on;
plot(tb,X1(5,Ta:T-1));hold on;
q=xlabel('${\rm Time (s)}$','interpreter','latex');
set(q,'Interpreter','latex');
q=ylabel('$$x(t)$$');
set(q,'Interpreter','latex');
q=legend('$$\Delta f(t)$$','$$\Delta X_{\rm t}(t)$$','$$\Delta X_{\rm g}(t)$$','$$\Delta X_{\rm e}(t)$$','$$\int {\rm ACE}(t)$$');
set(q,'Interpreter','latex');


