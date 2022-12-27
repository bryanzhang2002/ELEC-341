% Solution to ELEC 341 Assignment 2
% NAME: Bryan Zhang
% SN: 69238335
% DATE: 2022-09-23
% EMAIL: bryan.zhang@alumni.ubc.ca

clear all; clc; close all;

SN = 69238335;
s = tf('s');
a2DSPlot(SN);
fig = gcf;
axObjs = fig.Children;
dataObjs = axObjs.Children;
xData = dataObjs(1).XData;
yData = dataObjs(1).YData;
plot(xData*10^(-3), yData, 'k'); grid on; hold on;

%% Q1
Q1.Tr = 0.138; % Rise Time (ms)
Q1.Tp = 0.381; % Peak Time (ms)
Q1.Ts = 0.756; % Settle time (ms)

FV = 30.85;
overshoot = 45.7694 - FV; 
Q1.Pos = overshoot/FV * 100; % Percent Overshoot (%)

%% Q2: Approximate with rise time
Q2.Z = sqrt(log(overshoot/FV)^2/(pi^2+log(overshoot/FV)^2));
beta = sqrt(1-Q2.Z^2);
Q2.Wn = (1/(beta*Q1.Tr*1e-3)) * (pi - atan(beta/Q2.Z));

Q2.G = sysID(FV, Q2.Z, Q2.Wn);
step(Q2.G, 1.5e-3);

%% Q3
beta = sqrt(1-Q2.Z^2);
Q3.Wn = pi/(Q1.Tp*1e-3 * beta);

Q3.G = sysID(FV, Q2.Z, Q3.Wn);
step(Q3.G, 1.5e-3);

%% Q4
q4_overshoot = 2/3 * overshoot;
Q4.Z = sqrt(log(q4_overshoot/FV)^2/(pi^2+log(q4_overshoot/FV)^2));
Q4.Wn = pi/(Q1.Tp*1e-3 * beta);

Q4.G = sysID(FV, Q4.Z, Q4.Wn);
step(Q4.G, 1.5e-3);

%% Q5
q5_overshoot = 1/3 * overshoot;
Q5.Z = sqrt(log(q5_overshoot/FV)^2/(pi^2+log(q5_overshoot/FV)^2));
Q5.Wn = pi/(Q1.Tp*1e-3 * beta);

Q5.G = sysID(FV, Q5.Z, Q5.Wn);
step(Q5.G, 1.5e-3);

%% Q6
syms s w t
Q6.Z = 1;
step_resp_s = FV*(w^2/(s*(s^2 + 2*w*s + w^2)));
step_resp_t = ilaplace(step_resp_s, t);
Q6.Wn = double(vpasolve(subs(step_resp_t, t, Q1.Ts*1e-3) == 0.98*FV, w,100));

Q6.G = sysID(FV, Q6.Z, Q6.Wn);
step(Q6.G, 1.5e-3);

% num = [1];
% den = [1 1];
% filt = tf(num, den);
% [yf1, tf1] = step(filt); plot(tf1, yf1); hold on;
% [yf2, tf2] = step(filt*filt); plot(tf2, yf2);
% [yf3, tf3] = step(filt*filt*filt); plot(tf3, yf3);

%% Q7
Q7.Tr = (Q1.Tr*1e-3+0.217367e-3)/2;   % 0.217367ms approximated from plot of Q3
Q7.Te = Q7.Tr - Q1.Tr*1e-3;
wn = (1/(beta*Q7.Tr)) * (pi - atan(beta/Q2.Z));

Q7.G = sysID(FV, Q2.Z, wn);
step(Q7.G, 1.5e-3);

xline(Q1.Tr*1e-3, '-', 'Tr'); xline(Q1.Tp*1e-3, '-', 'Tp'); xline(Q1.Ts*1e-3, '-', 'Ts'); 
% yline(FV, '-', 'FV'); 
yline(FV*0.98, '-.'); yline(FV*1.02, '-.', 'FV');
legend('raw data','q2','q3','q4','q5','q6','q7', '','','');
grid on


function [f] = sysID(FV, Z, Wn)
    num = [FV*Wn^2];
    den = [1 2*Z*Wn Wn^2];
    f = tf(num,den);
end