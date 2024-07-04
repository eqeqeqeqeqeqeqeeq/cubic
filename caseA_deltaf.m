clc
clear
close all
tic
%% 系统矩阵设置
n=4; %绯荤粺缁存暟
D=1; M=10; Tt=0.3; Tg=0.1; R=0.05; bbeta=21;   alpha1=1; 
      
Kp=0.05;    
Ki=0.05;
A=[-D/M,   1/M,     0,   0;
    0, -1/Tt,  1/Tt,   0;
    -1/(R*Tg),     0, -1/Tg,   0;
    bbeta,     0,     0,   0];
Ad=[0,   0,     0,    0;
    0,   0,     0,     0;
    -(alpha1*Kp*bbeta)/Tg,   0,  0,   -(alpha1*Ki)/Tg;
      0,   0,     0,     0];  %为Bi


%% 
tau=@(t)150*sin(t)*sin(t);
d=0.01; % 时间间隔
t_begin=0; 
t_end=20000; 
ta=-2000:d:t_begin; Ta=length(ta);%Ta=10/0.02+1=501
tb=d:d:t_end;%[0.02-0.02-15]750     
t=[ta,tb];%[-10-0.02-15] 1250+1
T=length(t);%1251

X1=zeros(n,T);% 创造全零数组zeros,zeros(n)为返回一个 n×n 的全零矩阵后面往里面添加数据
%  X2=zeros(n,T);
% X3=zeros(n,T);
% X4=zeros(n,T);%0?
X0=[-0.2; 0.1;-0.3; -0.2]; % 初始值
%% 
for i=1:Ta
    X1(:,i)=X0;
    % X2(:,i)=X0;
    % X3(:,i)=X0;
    %X4(:,i)=X0;
    % phi=0; %扰动
end
for p=Ta:T%1251
   % phi=0.2*sin(p-Ta);%p=Ta:T， p-Ta相当于0时刻
    Tau=fix((tau((p-Ta)*d))/d); %fix为向左取整,tau((p-Ta)为具体时间的时滞
    X1(:,p+1)=X1(:,p)+d*A*X1(:,p)+d*Ad*X1(:,p-Tau); 
end


figure;
plot(tb,X1(1,Ta:T-1));hold on; %hold on是当前轴及图像保持而不被刷新，准备接受此后将绘制的图形，多图共存
q=legend('$$\Delta f(t)$$');%在坐标区上添加图例
% q=legend('$$x_1(t)$$','$$x_2(t)$$','$$x_3(t)$$','$$x_4(t)$$');
set(q,'Interpreter','latex');%添加一个小方框标加标签
