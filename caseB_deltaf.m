clc
clear
close all
tic

n=5; 
D=1; M=10; Tt=0.3; Tg=0.1; R=0.05; bbeta=300; rhoe=1/21; kappae=1.0; Te=0.1; alpha1=0.5; alpha2=0.5;
      
     
Kp11= 0.05 ;    
Ki11=0.1;
Kp21=0.05; 
Ki21=0.1;

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
Ad11=[0,   0,     0,   0,  0;
    0,   0,     0,   0,  0;
    -(alpha1*Kp11*bbeta)/Tg,   0,  0,  0, -(alpha1*Ki11)/Tg;
    -(alpha2*kappae*Kp11*bbeta)/Te,   0,  0,  0, -(alpha2*kappae*Ki11)/Te;
    0,   0,     0,   0,  0];
Ad21=[0,   0,     0,   0,  0;
    0,   0,     0,   0,  0;
    -(alpha1*Kp21*bbeta)/Tg,   0,  0,  0, -(alpha1*Ki21)/Tg;
    -(alpha2*kappae*Kp21*bbeta)/Te,   0,  0,  0, -(alpha2*kappae*Ki21)/Te;
    0,   0,     0,   0,  0];

%% 
tau1=@(t)0.5*sin(t)*sin(t)+10;

d=0.02;
t_begin=0;
t_end=100;
ta=-50:d:t_begin; Ta=length(ta);
tb=d:d:t_end;     t=[ta,tb];
T=length(t);

X1=zeros(n,T);

X0=[0.2; 0.1; -0.3; -0.5; 0.1];
%% 负荷扰动下系统状态轨迹图
for i=1:Ta
  
    X1(:,i)=X0;
 
end
for p=Ta:T
    Tau1=fix((tau1((p-Ta)*d))/d);
    X1(:,p+1)=X1(:,p)+ 0.5*(d*A1*X1(:,p)+d*Ad11*X1(:,p-Tau1))+ 0.5*(d*A2*X1(:,p)+d*Ad21*X1(:,p-Tau1));%%系统表达式
   
end
figure;

plot(tb,X1(1,Ta:T-1));hold on;
