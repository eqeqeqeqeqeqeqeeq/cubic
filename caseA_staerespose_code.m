clc
clear
close all
tic
%% 系统矩阵设置
n=4; %绯荤粺缁存暟
D=1; M=10; Tt=0.3; Tg=0.1; R=0.05; bbeta=21;   alpha1=1; 
      
Kp=0.1;    
Ki=0.15;
% Kp2=0.0987;    
% Ki2=0.1832;
% Kp3=0.0656;    
% Ki3=0.1232;
% Kp4=0.0441;    
% Ki4=0.0854;
% Kp2=0.1;1
% Ki2=0.05;
A=[-D/M,   1/M,     0,   0;
    0, -1/Tt,  1/Tt,   0;
    -1/(R*Tg),     0, -1/Tg,   0;
    bbeta,     0,     0,   0];
Ad=[0,   0,     0,    0;
    0,   0,     0,     0;
    -(alpha1*Kp*bbeta)/Tg,   0,  0,   -(alpha1*Ki)/Tg;
      0,   0,     0,     0];  %为Bi
% Ad2=[0,   0,     0,    0;
%     0,   0,     0,     0;
%     -(alpha1*Kp2*bbeta)/Tg,   0,  0,   -(alpha1*Ki2)/Tg;
%       0,   0,     0,     0];  %为Bi
% Ad3=[0,   0,     0,    0;
%     0,   0,     0,     0;
%     -(alpha1*Kp3*bbeta)/Tg,   0,  0,   -(alpha1*Ki3)/Tg;
%       0,   0,     0,     0];  %为Bi
% Ad4=[0,   0,     0,    0;
%     0,   0,     0,     0;
%     -(alpha1*Kp4*bbeta)/Tg,   0,  0,   -(alpha1*Ki4)/Tg;
%       0,   0,     0,     0];  %为Bi
% F=[-1/M;0;0;0];
  
% A1=[-D/M,   1/M,     0,   1/M, 0;
%     0, -1/Tt,  1/Tt,   0,  0;
%     -1/(R*Tg),     0, -1/Tg,  0,  0;
%     -(rhoe*kappae)/Te,  0,  0, -1/Te, 0;
%     bbeta,     0,     0,   0,   0];
% A2=[-D/M,   1/M,     0,   1/M, 0;
%     0, -1/Tt,  eta/Tt,   0,  0;
%     -1/(R*Tg),     0, -eta/Tg,  0,  0;
%     -(rhoe*kappae)/Te,  0,  0, -1/Te, 0;
%     bbeta,     0,     0,   0,   0];
% Ad1=[0,   0,     0,   0,  0;
%     0,   0,     0,   0,  0;
%     -(alpha1*Kp1*bbeta)/Tg,   0,  0,  0, -(alpha1*Ki1)/Tg;
%     -(alpha2*kappae*Kp1*bbeta)/Te,   0,  0,  0, -(alpha2*kappae*Ki1)/Te;
%     0,   0,     0,   0,  0];
% Ad2=[0,   0,     0,   0,  0;
%     0,   0,     0,   0,  0;
%     -(alpha1*Kp2*bbeta)/Tg,   0,  0,  0, -(alpha1*Ki2)/Tg;
%     -(alpha2*kappae*Kp2*bbeta)/Te,   0,  0,  0, -(alpha2*kappae*Ki2)/Te;
%     0,   0,     0,   0,  0];
% F=[-1/M,0,0,0,0;
%     0,  0,0,0,0;
%     0,  0,0,0,0;
%     0,  0,0,0,0;
%     0,  0,0,0,0];
% C=[bbeta,0,0,0,0;
%     0,0,0,0,1];

%% 
tau=@(t)6.669*sin(t)*sin(t)+2;
% tau2=@(t)0.5*sin(t)*sin(t)+2.5;
% tau3=@(t)0.5*sin(t)*sin(t)+3.5;
% tau4=@(t)0.5*sin(t)*sin(t)+4.5;  % @(t)为函数（t）引用  时滞
d=0.01; % 时间间隔
t_begin=0; 
t_end=80; 
ta=-20:d:t_begin; Ta=length(ta);%Ta=10/0.02+1=501
tb=d:d:t_end;%[0.02-0.02-15]750     
t=[ta,tb];%[-10-0.02-15] 1250+1
T=length(t);%1251

X1=zeros(n,T);% 创造全零数组zeros,zeros(n)为返回一个 n×n 的全零矩阵后面往里面添加数据
%  X2=zeros(n,T);
% X3=zeros(n,T);
% X4=zeros(n,T);%0?
X0=[0.2; 0.3;-0.3; -0.15]; % 初始值
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
% for p=Ta:T%1251
%    % phi=0.2*sin(p-Ta);%p=Ta:T， p-Ta相当于0时刻
%     Tau=fix((tau((p-Ta)*d))/d); %fix为向左取整,tau((p-Ta)为具体时间的时滞
%     X1(:,p+1)=X2(:,p)+d*A*X2(:,p)+d*Ad*X2(:,p-Tau); 
% end
% for p=Ta:T%1251
%    % phi=0.2*sin(p-Ta);%p=Ta:T， p-Ta相当于0时刻
%     Tau=fix((tau((p-Ta)*d))/d); %fix为向左取整,tau((p-Ta)为具体时间的时滞
%     X3(:,p+1)=X3(:,p)+d*A*X3(:,p)+d*Ad*X3(:,p-Tau); 
% end
% for p=Ta:T%1251
%    % phi=0.2*sin(p-Ta);%p=Ta:T， p-Ta相当于0时刻
%     Tau=fix((tau((p-Ta)*d))/d); %fix为向左取整,tau((p-Ta)为具体时间的时滞
%     X4(:,p+1)=X4(:,p)+d*A*X4(:,p)+d*Ad*X4(:,p-Tau); 
% end

figure;
% subplot(3,1,1) %subplot(m,n,p) 将当前图窗划分为 m×n 网格，并在 p 指定的位置创建坐标区。
plot(tb,X1(1,Ta:T-1));hold on; %hold on是当前轴及图像保持而不被刷新，准备接受此后将绘制的图形，多图共存
plot(tb,X1(2,Ta:T-1));hold on; %画下一个图时，保持本次图
plot(tb,X1(3,Ta:T-1));hold on;%1 2 3 4 WHAT?
 plot(tb,X1(4,Ta:T-1));hold on;%tb-xlabel X(4,Ta:T-1)-ylabel 矩阵元素一一对应绘制第四行即第四个变量的曲线
% plot(tb,X(5,Ta:T-1));hold on;
q=xlabel('${\rm Time (s)}$','interpreter','latex');%为x轴添加标签
set(q,'Interpreter','latex'); 
q=ylabel('$$x(t)$$');
 set(q,'Interpreter','latex');%为Y轴添加标签
q=legend('$$\Delta f(t)$$','$$\Delta Z_{t}(t)$$','$$\Delta Z_{g}(t)$$','$$\int {\rm ACE(t)}$$');%在坐标区上添加图例
% q=legend('$$x_1(t)$$','$$x_2(t)$$','$$x_3(t)$$','$$x_4(t)$$');
set(q,'Interpreter','latex');%添加一个小方框标加标签
% h=legend('$$\sqrt{x^2+y^2}$$');
% set(h,'Interpreter','latex')在legend里加数学字符$数学公式$ 或者 $$数学公式$$

% subplot(3,1,2)
% figure;
% plot(tb,X2(1,Ta:T-1));hold on;
% plot(tb,X2(2,Ta:T-1));hold on;
% plot(tb,X2(3,Ta:T-1));hold on;
% plot(tb,X2(4,Ta:T-1));hold on;
% q=xlabel('${\rm Time (s)}$','interpreter','latex');
% set(q,'Interpreter','latex');
% q=ylabel('$$x(t)$$');
% set(q,'Interpreter','latex');
% q=legend('$$\Delta f(t)$$','$$\Delta P_{\rm t}(t)$$','$$\Delta X_{\rm G}(t)$$','$$\int {\rm ACE}(t)$$');
% % q=legend('$$x_1(t)$$','$$x_2(t)$$','$$x_3(t)$$','$$x_4(t)$$');
% set(q,'Interpreter','latex');
% % subplot(3,1,3)
 %figure;
%plot(tb,X3(1,Ta:T-1));hold on;
%plot(tb,X3(2,Ta:T-1));hold on;
%plot(tb,X3(3,Ta:T-1));hold on;
%plot(tb,X3(4,Ta:T-1));hold on;
%q=xlabel('${\rm Time (s)}$','interpreter','latex');
%set(q,'Interpreter','latex');
%q=ylabel('$$x(t)$$');
%set(q,'Interpreter','latex');
%q=legend('$$\Delta f(t)$$','$$\Delta P_{\rm t}(t)$$','$$\Delta X_{\rm G}(t)$$','$$\int {\rm ACE}(t)$$');
% q=legend('$$x_1(t)$$','$$x_2(t)$$','$$x_3(t)$$','$$x_4(t)$$');
%set(q,'Interpreter','latex');
% subplot(2,2,4)
% % figure;
% plot(tb,X4(1,Ta:T-1));hold on;
% plot(tb,X4(2,Ta:T-1));hold on;
% plot(tb,X4(3,Ta:T-1));hold on;
% plot(tb,X4(4,Ta:T-1));hold on;
% q=xlabel('${\rm Time (s)}$','interpreter','latex');
% set(q,'Interpreter','latex');
% q=ylabel('$$x(t)$$');
% set(q,'Interpreter','latex');
% q=legend('$$\Delta f(t)$$','$$\Delta P_{\rm t}(t)$$','$$\Delta X_{\rm G}(t)$$','$$\int {\rm ACE}(t)$$');
% % q=legend('$$x_1(t)$$','$$x_2(t)$$','$$x_3(t)$$','$$x_4(t)$$');
% set(q,'Interpreter','latex');
