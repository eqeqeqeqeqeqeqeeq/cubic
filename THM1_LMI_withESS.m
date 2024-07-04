%-------------------------------------------------------------------
%--------     cublic function  nagtive defination     ----------------
clc
clear
close all
tic
%Kp=0.1，Ki=0.05 11.989 ;Kp=0.1，Ki=0.1 10.795; Kp=0.1，Ki=0.15 8.671;
%% 系统矩阵设置

n=5; %系统维数
D=1; M=10; Tt=0.3; Tg=0.1; R=0.05; bbeta=21; rhoe=1/21; kappae=1.0; Te=0.1; alpha1=0.5; alpha2=0.5;
miu1=0.5;
miu2=0.5;
     
Kp11= 0.10;    
Ki11=0.15;  
Kp21=0.10;
Ki21=0.15;

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
A=(miu1)*A1+(miu2)*A2;%dl=1,zt=0.2
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
Ad=0.5*Ad11+0.5*Ad21;
h1=0;ACC1=h1;
ACC2 = 25; acc = 0.001;
%Kp=0.1  ;  Ki=0.1; 12.692
%Kp=0.1  ;  Ki=0.15;%9.210
while true
h2 = (ACC1 + ACC2) / 2;
h12 = h2 - h1; %tao2-tao1
%% LMI变量定义
num=12;
for i=1:num
    e=[zeros(n,(i-1)*n), eye(n), zeros(n,(num-i)*n)]; % 定义ei
    eval(['e',num2str(i),'=[e]']); %计算 MATLAB 表达式，即计算e
end
es=A*e1+Ad*e3; 
setlmis([])
[P  P_n   sP]=lmivar(1,[7*n,1]);
[Q1 Q1_n sQ1]=lmivar(1,[2*n,1]);
[Q2 Q2_n sQ2]=lmivar(1,[6*n,1]);
[R1 R1_n sR1]=lmivar(1,[n,1]);
[R2 R2_n sR2]=lmivar(1,[n,1]);
[N1 N1_n sN1]=lmivar(1,[3*n,1]);
[N5 N5_n sN5]=lmivar(1,[3*n,1]);
[N2 N2_n sN2]=lmivar(2,[3*n,3*n]);
[N4 N4_n sN4]=lmivar(2,[3*n,3*n]);
[N3 N3_n sN3]=lmivar(2,[3*n,5*n]);
hatR1=lmivar(3,[sR1,        zeros(n,n),  zeros(n,n);
                zeros(n,n), sR1,         zeros(n,n);
                zeros(n,n), zeros(n,n),       sR1]);
hatR2=lmivar(3,[sR2,        zeros(n,n),  zeros(n,n);
                zeros(n,n), sR2,         zeros(n,n);
                zeros(n,n), zeros(n,n),  sR2]);
L1=[e1;e2;e4;h1*e5;-h1*e6+h2*e7;h1*h1*e8;h1*h1*e9+h2*h2*e10-h1*h2*e6];
L2=[zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);e6-e7;zeros(n,num*n);-2*h1*e9-2*h2*e10+(h1+h2)*e6];
L3=[zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);e9+e10-e6];
L4=[es;e11;e12;e1-e2;e2-e4;h1*(e1-e5);h12*e2+h1*e6-h2*e7];
L5=[zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);-e6+e7];
L6=[e2;e11;e1;e2;e4;zeros(n,num*n)];
L7=[e4;e12;e1;e2;e4;-h1*e6+h2*e7];
L8=[-h1*e6+h2*e7;e2-e4;h12*e1;h12*e2;h12*e4;h1*h1*e9+h2*h2*e10-h1*h2*e6];
L9=[e6-e7;zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);-2*h1*e9-2*h2*e10+(h1+h2)*e6];
L10=[zeros(n,num*n);zeros(n,num*n);es;e11;e12;e2];
L11=[zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);e9+e10-e6];
L12=[zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);-e6+e7];
gamma0=[e7;e8;e9;e11;e12];
gamma1=[e1-e2; e1+e2-2*e5; e1-e2+6*e5-12*e8];
gamma2=[e2-e3; e2+e3-2*e6; e2-e3+6*e6-12*e9];
gamma3=[e3-e4; e3+e4-2*e7; e3-e4+6*e7-12*e10];
Coe=[eye(n),      zeros(n),         zeros(n);
    zeros(n),    sqrt(3)*eye(n),   zeros(n);
    zeros(n),    zeros(n),         sqrt(5)*eye(n)];
%% f(h1)
i=1;
%f(h1)
lmiterm([i,1,1,P],L1'+h1*L2'+h1*h1*L3',L4+h1*L5,'s');
lmiterm([i,1,1,Q1],[e1',es'],[e1',es']');
lmiterm([i,1,1,Q1],-[e2',e11'],[e2',e11']');
lmiterm([i,1,1,Q2],L6',L6);
lmiterm([i,1,1,Q2],-L7'+h1*L12',L7-h1*L12);
lmiterm([i,1,1,R1],es'*h1*h1,es);
lmiterm([i,1,1,R2],es'*h12*h12,es);
lmiterm([i,1,1,Q2],L8'+h1*L9'+h1*h1*L11',L10,'s');
lmiterm([i,1,1,hatR1],-gamma1'*Coe,Coe*gamma1);
lmiterm([i,1,1,hatR2],-2*gamma2'*Coe,Coe*gamma2);
lmiterm([i,1,1,hatR2],-gamma3'*Coe,Coe*gamma3);
lmiterm([i,1,1,N1],-gamma2',gamma2,'s');
lmiterm([i,1,1,N2],-gamma2',gamma3,'s');
lmiterm([i,1,1,N3],-gamma2',gamma0,'s');
%Ngamma3'
lmiterm([i,1,2,-N4],gamma2',1);
lmiterm([i,1,2,-N5],gamma3',1);
lmiterm([i,1,2,-N3],-gamma0',1);

lmiterm([i,2,2,hatR2],-Coe,Coe);
%% f(h2)
i=2;
%fai(h2)
lmiterm([i,1,1,P],L1'+h2*L2'+h2*h2*L3',L4+h2*L5,'s');
lmiterm([i,1,1,Q1],[e1',es'],[e1',es']');
lmiterm([i,1,1,Q1],-[e2',e11'],[e2',e11']');
lmiterm([i,1,1,Q2],L6',L6);
lmiterm([i,1,1,Q2],-L7'+h2*L12',L7-h2*L12);
lmiterm([i,1,1,R1],es'*h1*h1,es);
lmiterm([i,1,1,R2],es'*h12*h12,es);
lmiterm([i,1,1,Q2],L8'+h2*L9'+h2*h2*L11',L10,'s');
lmiterm([i,1,1,hatR1],-gamma1'*Coe,Coe*gamma1);
lmiterm([i,1,1,hatR2],-gamma2'*Coe,Coe*gamma2);
lmiterm([i,1,1,hatR2],-2*gamma3'*Coe,Coe*gamma3);
lmiterm([i,1,1,N4],-gamma3',gamma2,'s');
lmiterm([i,1,1,N5],-gamma3',gamma3,'s');
lmiterm([i,1,1,N3],gamma3',gamma0,'s');
%Ngamma2'
lmiterm([i,1,2,-N1],gamma2',1);
lmiterm([i,1,2,-N2],gamma3',1);
lmiterm([i,1,2,-N3],gamma0',1);

lmiterm([i,2,2,hatR2],-Coe,Coe);
%% f(h1)-a2*h1^2
i=3;
%f(h1)
lmiterm([i,1,1,P],L1'+h1*L2'+h1*h1*L3',L4+h1*L5,'s');
lmiterm([i,1,1,Q1],[e1',es'],[e1',es']');
lmiterm([i,1,1,Q1],-[e2',e11'],[e2',e11']');
lmiterm([i,1,1,Q2],L6',L6);
lmiterm([i,1,1,Q2],-L7'+h1*L12',L7-h1*L12);
lmiterm([i,1,1,R1],es'*h1*h1,es);
lmiterm([i,1,1,R2],es'*h12*h12,es);
lmiterm([i,1,1,Q2],L8'+h1*L9'+h1*h1*L11',L10,'s');
lmiterm([i,1,1,hatR1],-gamma1'*Coe,Coe*gamma1);
lmiterm([i,1,1,hatR2],-2*gamma2'*Coe,Coe*gamma2);
lmiterm([i,1,1,hatR2],-gamma3'*Coe,Coe*gamma3);
lmiterm([i,1,1,N1],-gamma2',gamma2,'s');
lmiterm([i,1,1,N2],-gamma2',gamma3,'s');
lmiterm([i,1,1,N3],-gamma2',gamma0,'s');
%-a2*h1^2
lmiterm([i,1,1,P],-h1*h1*L2',L5,'s');
lmiterm([i,1,1,P],-h1*h1*L3',L4,'s');
lmiterm([i,1,1,Q2],-h1*h1*L11',L10,'s');
lmiterm([i,1,1,Q2],h1*h1*L12',L12);

%Ngamma3'
lmiterm([i,1,2,-N4],gamma2',1);
lmiterm([i,1,2,-N5],gamma3',1);
lmiterm([i,1,2,-N3],-gamma0',1);

lmiterm([i,2,2,hatR2],-Coe,Coe);
%% f(h2)-a2*h2^2
i=4;
%fai(h2)
lmiterm([i,1,1,P],L1'+h2*L2'+h2*h2*L3',L4+h2*L5,'s');
lmiterm([i,1,1,Q1],[e1',es'],[e1',es']');
lmiterm([i,1,1,Q1],-[e2',e11'],[e2',e11']');
lmiterm([i,1,1,Q2],L6',L6);
lmiterm([i,1,1,Q2],-L7'+h2*L12',L7-h2*L12);
lmiterm([i,1,1,R1],es'*h1*h1,es);
lmiterm([i,1,1,R2],es'*h12*h12,es);
lmiterm([i,1,1,Q2],L8'+h2*L9'+h2*h2*L11',L10,'s');
lmiterm([i,1,1,hatR1],-gamma1'*Coe,Coe*gamma1);
lmiterm([i,1,1,hatR2],-gamma2'*Coe,Coe*gamma2);
lmiterm([i,1,1,hatR2],-2*gamma3'*Coe,Coe*gamma3);
lmiterm([i,1,1,N4],-gamma3',gamma2,'s');
lmiterm([i,1,1,N5],-gamma3',gamma3,'s');
lmiterm([i,1,1,N3],gamma3',gamma0,'s');
%-a2*h2^2
lmiterm([i,1,1,P],-h2*h2*L2',L5,'s');
lmiterm([i,1,1,P],-h2*h2*L3',L4,'s');
lmiterm([i,1,1,Q2],-h2*h2*L11',L10,'s');
lmiterm([i,1,1,Q2],h2*h2*L12',L12);
%Ngamma2^T
lmiterm([i,1,2,-N1],gamma2',1);
lmiterm([i,1,2,-N2],gamma3',1);
lmiterm([i,1,2,-N3],gamma0',1);

lmiterm([i,2,2,hatR2],-Coe,Coe);
%% fai(h1)-h1^3*a3
i=5;
%f(h1)
lmiterm([i,1,1,P],L1'+h1*L2'+h1*h1*L3',L4+h1*L5,'s');
lmiterm([i,1,1,Q1],[e1',es'],[e1',es']');
lmiterm([i,1,1,Q1],-[e2',e11'],[e2',e11']');
lmiterm([i,1,1,Q2],L6',L6);
lmiterm([i,1,1,Q2],-L7'+h1*L12',L7-h1*L12);
lmiterm([i,1,1,R1],es'*h1*h1,es);
lmiterm([i,1,1,R2],es'*h12*h12,es);
lmiterm([i,1,1,Q2],L8'+h1*L9'+h1*h1*L11',L10,'s');
lmiterm([i,1,1,hatR1],-gamma1'*Coe,Coe*gamma1);
lmiterm([i,1,1,hatR2],-2*gamma2'*Coe,Coe*gamma2);
lmiterm([i,1,1,hatR2],-gamma3'*Coe,Coe*gamma3);
lmiterm([i,1,1,N1],-gamma2',gamma2,'s');
lmiterm([i,1,1,N2],-gamma2',gamma3,'s');
lmiterm([i,1,1,N3],-gamma2',gamma0,'s');
%-h1^3*a3
lmiterm([i,1,1,P],-h1*h1*h1*L3',L5,'s');

%Ngamma3'
lmiterm([i,1,2,-N4],gamma2',1);
lmiterm([i,1,2,-N5],gamma3',1);
lmiterm([i,1,2,-N3],-gamma0',1);

lmiterm([i,2,2,hatR2],-Coe,Coe);

%% fai(h2)-h2^3*a3
i=6;
%fai(h2)
lmiterm([i,1,1,P],L1'+h2*L2'+h2*h2*L3',L4+h2*L5,'s');
lmiterm([i,1,1,Q1],[e1',es'],[e1',es']');
lmiterm([i,1,1,Q1],-[e2',e11'],[e2',e11']');
lmiterm([i,1,1,Q2],L6',L6);
lmiterm([i,1,1,Q2],-L7'+h2*L12',L7-h2*L12);
lmiterm([i,1,1,R1],es'*h1*h1,es);
lmiterm([i,1,1,R2],es'*h12*h12,es);
lmiterm([i,1,1,Q2],L8'+h2*L9'+h2*h2*L11',L10,'s');
lmiterm([i,1,1,hatR1],-gamma1'*Coe,Coe*gamma1);
lmiterm([i,1,1,hatR2],-gamma2'*Coe,Coe*gamma2);
lmiterm([i,1,1,hatR2],-2*gamma3'*Coe,Coe*gamma3);
lmiterm([i,1,1,N4],-gamma3',gamma2,'s');
lmiterm([i,1,1,N5],-gamma3',gamma3,'s');
lmiterm([i,1,1,N3],gamma3',gamma0,'s');
%-h2^3*a3
lmiterm([i,1,1,P],-h2*h2*h2*L3',L5,'s');

%Ngamma2^T
lmiterm([i,1,2,-N1],gamma2',1);
lmiterm([i,1,2,-N2],gamma3',1);
lmiterm([i,1,2,-N3],gamma0',1);

lmiterm([i,2,2,hatR2],-Coe,Coe);

%% a1*h1+a0
i=7;
%a1'*h1
lmiterm([i,1,1,P],h1*L2',L4,'s');
lmiterm([i,1,1,P],h1*L1',L5,'s');
lmiterm([i,1,1,Q2],h1*L9',L10,'s');
lmiterm([i,1,1,Q2],h1*L7',L12,'s');
lmiterm([i,1,1,N1],h1*(1/h12)*gamma2',gamma2,'s');
lmiterm([i,1,1,N2],h1*(1/h12)*gamma2',gamma3,'s');
lmiterm([i,1,1,N3],h1*(1/h12)*gamma2',gamma0,'s');
lmiterm([i,1,1,N4],h1*(-1/h12)*gamma3',gamma2,'s');
lmiterm([i,1,1,N5],h1*(-1/h12)*gamma3',gamma3,'s');
lmiterm([i,1,1,N3],h1*(-1/h12)*gamma3'*(-1),gamma0,'s');
lmiterm([i,1,1,hatR2],h1*(1/h12)*gamma2'*Coe,Coe*gamma2);
lmiterm([i,1,1,hatR2],h1*(-1/h12)*gamma3'*Coe,Coe*gamma3);
%a0'
lmiterm([i,1,1,P],L1',L4,'s');
lmiterm([i,1,1,Q2],L8',L10,'s');
lmiterm([i,1,1,N1],(-h2/h12)*gamma2',gamma2,'s');
lmiterm([i,1,1,N2],(-h2/h12)*gamma2',gamma3,'s');
lmiterm([i,1,1,N3],(-h2/h12)*gamma2',gamma0,'s');
lmiterm([i,1,1,N4],(h1/h12)*gamma3',gamma2,'s');
lmiterm([i,1,1,N5],(h1/h12)*gamma3',gamma3,'s');
lmiterm([i,1,1,N3],(h1/h12)*gamma3'*(-1),gamma0,'s');
lmiterm([i,1,1,Q1],[e1',es'],[e1',es']');
lmiterm([i,1,1,Q1],-[e2',e11'],[e2',e11']');
lmiterm([i,1,1,Q2],L6',L6);
lmiterm([i,1,1,Q2],-L7',L7);
lmiterm([i,1,1,R1],es'*h1*h1,es);
lmiterm([i,1,1,R2],es'*h12*h12,es);
lmiterm([i,1,1,hatR1],-gamma1'*Coe,Coe*gamma1);
lmiterm([i,1,1,hatR2],-(2+(h1/h12))*gamma2'*Coe,Coe*gamma2);
lmiterm([i,1,1,hatR2],-(1-(h1/h12))*gamma3'*Coe,Coe*gamma3);

%Ngamma3^T
lmiterm([i,1,2,-N4],gamma2',1);
lmiterm([i,1,2,-N5],gamma3',1);
lmiterm([i,1,2,-N3],-gamma0',1);

lmiterm([i,2,2,hatR2],-Coe,Coe);%t>=2

%% a1*h2+a0
i=8;
%a1'*h2
lmiterm([i,1,1,P],h2*L2',L4,'s');
lmiterm([i,1,1,P],h2*L1',L5,'s');
lmiterm([i,1,1,Q2],h2*L9',L10,'s');
lmiterm([i,1,1,Q2],h2*L7',L12,'s');
lmiterm([i,1,1,N1],h2*(1/h12)*gamma2',gamma2,'s');
lmiterm([i,1,1,N2],h2*(1/h12)*gamma2',gamma3,'s');
lmiterm([i,1,1,N3],h2*(1/h12)*gamma2',gamma0,'s');
lmiterm([i,1,1,N4],h2*(-1/h12)*gamma3',gamma2,'s');
lmiterm([i,1,1,N5],h2*(-1/h12)*gamma3',gamma3,'s');
lmiterm([i,1,1,N3],h2*(-1/h12)*gamma3'*(-1),gamma0,'s');
lmiterm([i,1,1,hatR2],h2*(1/h12)*gamma2'*Coe,Coe*gamma2);
lmiterm([i,1,1,hatR2],h2*(-1/h12)*gamma3'*Coe,Coe*gamma3);
%a0'
lmiterm([i,1,1,P],L1',L4,'s');
lmiterm([i,1,1,Q2],L8',L10,'s');
lmiterm([i,1,1,N1],(-h2/h12)*gamma2',gamma2,'s');
lmiterm([i,1,1,N2],(-h2/h12)*gamma2',gamma3,'s');
lmiterm([i,1,1,N3],(-h2/h12)*gamma2',gamma0,'s');
lmiterm([i,1,1,N4],(h1/h12)*gamma3',gamma2,'s');
lmiterm([i,1,1,N5],(h1/h12)*gamma3',gamma3,'s');
lmiterm([i,1,1,N3],(h1/h12)*gamma3'*(-1),gamma0,'s');
lmiterm([i,1,1,Q1],[e1',es'],[e1',es']');
lmiterm([i,1,1,Q1],-[e2',e11'],[e2',e11']');
lmiterm([i,1,1,Q2],L6',L6);
lmiterm([i,1,1,Q2],-L7',L7);
lmiterm([i,1,1,R1],es'*h1*h1,es);
lmiterm([i,1,1,R2],es'*h12*h12,es);
lmiterm([i,1,1,hatR1],-gamma1'*Coe,Coe*gamma1);
lmiterm([i,1,1,hatR2],-(2+(h1/h12))*gamma2'*Coe,Coe*gamma2);
lmiterm([i,1,1,hatR2],-(1-(h1/h12))*gamma3'*Coe,Coe*gamma3);

%Ngamma2^T
lmiterm([i,1,2,-N1],gamma2',1);
lmiterm([i,1,2,-N2],gamma3',1);
lmiterm([i,1,2,-N3],gamma0',1);


lmiterm([i,2,2,hatR2],-Coe,Coe);%t>=2

%%
i=i+1;

lmiterm([-i,1,1,P],1,1);

i=i+1;
lmiterm([-i,1,1,Q1],1,1);

i=i+1;
lmiterm([-i,1,1,Q2],1,1);

i=i+1;
lmiterm([-i,1,1,R1],1,1);

i=i+1;
lmiterm([-i,1,1,R2],1,1);
%% 求解语句

lmisys=getlmis;
[tmin,xfeas]=feasp(lmisys);
if tmin < 0
    ACC1 = h2;
    if ACC2 - ACC1 <= acc
        break;
    end
else 
    ACC2 = h2;

 end

end
