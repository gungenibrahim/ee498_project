addpath(fullfile(matlabroot,'examples','mpc_featured','main'));
mdl = 'mpcLKAsystem';
open_system(mdl)
Ts = 0.1;
T = 15;
% System parameters
m = 2023;       %kg
lf = 1.26;      %m
Cf = 2.864e5;   %N/rad
Iz = 6286;       %kg.m2
lr = 1.90;      %m
Cr = 1.948e5;   %N/rad

%% changable values    
x0 = [ 0;0;0;10];
p_dis = 0.001; %1/R
v = 30;         %m/sec
ls = 10; % m

Vx = v
%%


a11 = -(Cf+Cr)/(m*v);
a12 = -1-(Cf*lf-Cr*lr)/(m*v^2);
a21 = -(Cf*lf-Cr*lr)/(Iz);
a22 = -(Cf*lf^2 + Cr*lr^2)/(Iz*v);
b1 = Cf / (m*v);
b2 = Cf*lf / Iz;



% Finite Horizon MPC as quadratic program
% State-Space Model of the System
A = [a11 a12 0 0;
      a21 a22 0 0; 
      0 1 0 0; 
      v ls v 0];
B = [b1; b2; 0; 0];
C =  [0 0 0 1];
D =  [0];


G = ss(A,B,C,D);

PredictionHorizon = 4;

time = 0:0.1:15;
md = getCurvature(Vx,time);

u_min = -0.349;
u_max = 0.349;

sim(mdl)
mpcLKAplot(logsout)

