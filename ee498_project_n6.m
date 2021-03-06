%% EE498 Term Project
%%
% �brahim G�ngen
% 1936939

%%
% Some figure formatting
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex'); 
clc
clear all

%% Introduction
%%
% 
%
%%


% System parameters
m = 2023;       %kg
lf = 1.26;      %m
cf = 2.864e5;   %N/rad
J = 6286;       %kg.m2
lr = 1.90;      %m
cr = 1.948e5;   %N/rad

%% changable values    
x0 = [ 0;0;0;10];
p_dis = 0.001; %1/R
v = 30;         %m/sec
ls = 10; % m
%%


a11 = -(cf+cr)/(m*v);
a12 = -1-(cf*lf-cr*lr)/(m*v^2);
a21 = -(cf*lf-cr*lr)/(J);
a22 = -(cf*lf^2 + cr*lr^2)/(J*v);
b1 = cf / (m*v);
b2 = cf*lf / J;



% Finite Horizon MPC as quadratic program
% State-Space Model of the System
Ac = [a11 a12 0 0;
      a21 a22 0 0; 
      0 1 0 0; 
      v ls v 0];
Bc = [b1; b2; 0; 0];
C =  [0 0 0 1];
D =  [0];

Ts = 0.05;

%%
T = sym('T');
x = sym('x');

A_dis = expm(Ac*T);
B_dis = int(expm(Ac*x),x,[0 T])*Bc;
%%

A_dis = subs(A_dis,T,Ts);
B_dis = subs(B_dis,T,Ts);
A = double(A_dis);
B = real(double(B_dis));


%  B_dis = (int(expm(Ac*x),x,[0,T])*Bc)
% % 
%  A = subs(A,T,Ts)
%  B = subs(B,T,Ts)
% 
% A = [   0.449,      -0.0492,     0,      0
%         0.0728,      0.537,      0,      0
%         0.00467,     0.0746,    1.0,     0
%         2.133,       0.998,     3.0,    1.0]
%     
% B =  [  0.144
%         4.303
%         0.237
%         3.677]
    


% Optimal control solution for $N = 4$
G = [zeros(4,1) zeros(4,1) zeros(4,1) zeros(4,1) zeros(4,1) zeros(4,1); 
    B zeros(4,1) zeros(4,1) zeros(4,1) zeros(4,1) zeros(4,1);
    A*B B zeros(4,1) zeros(4,1) zeros(4,1) zeros(4,1); 
    A^2*B A*B B zeros(4,1) zeros(4,1) zeros(4,1);
    A^3*B A^2*B A*B B zeros(4,1) zeros(4,1);
    A^4*B A^3*B A^2*B A*B B zeros(4,1)
    A^5*B A^4*B A^3*B A^2*B A*B B];
H = [eye(4); A; A^2; A^3; A^4; A^5; A^6];
Q = C'*C;
R = 0.001;

% Q = eye(2);
Pinf = dare(A,B,Q,R,zeros(4,1),eye(4));
Kinf = inv(R+B'*Pinf*B)*B'*Pinf*A;
% A*X*A' - X + Q = 0;  X = dlyap(A,Q)
P = dlyap( (A-B*Kinf)',Q+Kinf'*R*Kinf );
Qf = P;
Qbar = blkdiag(Q,Q,Q,Q,Q,Q,Qf);
Rbar = blkdiag(R,R,R,R,R,R);
M = G'*Qbar*G + Rbar;
% input bound: umin <= u <= umax
umin = -0.3491; % -20 degrees
umax = 0.3491; % 20 degrees
lb = [umin;umin;umin;umin];
ub = [umax;umax;umax;umax];
% Apply MPC steps
xVec(:,1) = x0;
yVec(1) = C*x0;
uVec = [];
for kk = 1:60
    alpha = G'*Qbar*H*xVec(:,kk);
    Usol = quadprog(M,alpha',[],[],[],[],lb,ub);
    uVec(kk) = Usol(1);
    xVec(:,kk+1) = A*xVec(:,kk) + B*uVec(kk) + [0; 0; -v*p_dis; 0];
    yVec(kk+1) = C*xVec(:,kk+1);
    Xsol(:,1) = xVec(:,kk);
    Xsol(:,2) = A*Xsol(:,1) + B*Usol(1);
    Xsol(:,3) = A*Xsol(:,2) + B*Usol(2);
    Xsol(:,4) = A*Xsol(:,3) + B*Usol(3);
    Ysol(1) = C*Xsol(:,1);
    Ysol(2) = C*Xsol(:,2);
    Ysol(3) = C*Xsol(:,3);
    Ysol(4) = C*Xsol(:,4);
end

uVec = [uVec uVec(end)];
tVec = [0:Ts:3];
%figure;
subplot(3,2,1)
stairs(tVec,uVec,'LineWidth',2);
hold all;
xlabel('time [sec]')
grid
ylabel('$u$')
title('Input $u$')

subplot(3,2,2)
stairs(tVec,C*xVec,'LineWidth',2)
hold all;
grid
xlabel('time [sec]')
ylabel('$y$')
title('Output $y$')

subplot(3,2,3)
stairs(tVec,[1 0 0 0]*xVec,'LineWidth',2)
hold all;
grid
xlabel('time [sec]')
ylabel('$x_1$')
title('State $x_1$')

subplot(3,2,4)
stairs(tVec,[0 1 0 0]*xVec,'LineWidth',2)
hold all;
grid
xlabel('time [sec]')
ylabel('$x_2$')
title('State $x_2$')

subplot(3,2,5)
stairs(tVec,[0 0 1 0]*xVec,'LineWidth',2)
hold all;
grid
xlabel('time [sec]')
ylabel('$x_3$')
title('State $x_3$')

subplot(3,2,6)
stairs(tVec,[0 0 0 1]*xVec,'LineWidth',2)
hold all;
grid
xlabel('time [sec]')
ylabel('$x_4$')
title('State $x_4$')

set(findall(gcf,'Type','line'),'LineWidth',2)
set(findall(gcf,'-property','FontSize'),'FontSize',14);
% legend('$u_{max} = 1.5$','$u_{max} = 2.5$','$u_{max} = 4$')

