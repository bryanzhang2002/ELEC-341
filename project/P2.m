% Solution to ELEC 341 Project 2
% NAME: Bryan Zhang
% SN: 69238335
% DATE: 2022-12-04
% EMAIL: bryan.zhang@alumni.ubc.ca
clear all; close all; clc; load part1.mat

SN = 69238335; sn = [6 9 2 3 8 3 3 5];
A = sn(1)+10; B = sn(2)+10; C = sn(3)+10; D = sn(4)+10; E = sn(5)+10; F = sn(6)+10; G = sn(7)+10; H = sn(8)+10;
s = tf('s');

%% Q8
Q8.N = 6;
Q8.Nhat = 1.94;
Q8.Dp = (1/s)*2*CF/(Q8.Nhat*s+2*CF);

%% Q9
Q9.Kref = margin(Q8.Dp*Q5.GH);
[~, ~, Q9.wxo] = margin(Q9.Kref*Q8.Dp*Q5.GH);

%% Q10
sigma = -11.7264; omega = 3.23;
z1 = sigma - omega*1j;
z2 = sigma + omega*1j;

Q10.Z = [z1 z2];
Q10.PM = 127.3465;
Q10.D = 1/(z1*z2)*(s-z1)*(s-z2)*Q8.Dp;

%% real part
pm_max = 0;
omega  = 3.23;
for sigma=Q9.wxo-2:0.01:Q9.wxo+2
    sigma = -sigma; disp(sigma)
    z1 = sigma-omega*1i;
    z2 = sigma+omega*1i;
    Q10D = 1/(z1*z2)*(s-z1)*(s-z2)*Q8.Dp;
    [~,pm,~,~] = margin(Q5.GH*Q10D*Q9.Kref);
    if(pm>=pm_max)
        pm_max = pm;
        value_sigma = sigma;
    end
end
pm_max
value_sigma

%% imaginary part
pm_max = 0;
sigma = -11.7264;
for omega=0:0.01:5
    disp(omega)
    z1 = sigma-omega*1i;
    z2 = sigma+omega*1i;
    Q10D = 1/(z1*z2)*(s-z1)*(s-z2)*Q8.Dp;
    [~,pm,~,~] = margin(Q5.GH*Q10D*Q9.Kref);
    if(pm>=pm_max)
        pm_max = pm;
        value_omega = omega;
    end
end
pm_max
value_omega

%% Q11
p = -2*CF/Q8.Nhat;
Q11.K = 3.7009;
Q11.Kp = 1/p - (z1+z2)/(z1*z2);
Q11.Ki = 1;
Q11.Kd = 1/(z1*z2) + Q11.Kp/p;

%%  manually tune the bounds
for Q11K = 3.7:0.0001:3.701
    disp(Q11K)
    [Gm,Pm,Wcg,Wcp] = margin(Q5.GH*Q10.D*Q11K);
    disp(Pm)
    if(Pm>=29.9 && Pm<=30.1)    % manually tune this
        Pm_30 = Pm;
        Q11.K = Q11K;
    end
end
Pm_30

%% Q12
%{
Requirments:
    1. OS < 20%
    2. Ts < 200ms
Constraints:
    1. you may not change the control frequency
    2. you may not change the derivative filter
Goals (in decreasing order of importance)
    1. No undershoot (never dips below FV after peak)
    2. OS as small as possible
    3. Ts as small as possible
%}
% fuck it just do this randomly no cap
close all
tune1 = 0.54;    % tune K
tune2 = 0.035;    % tune Kd
tune3 = 1.45;     % tune Kp
factor_for_K = 0.45;    % tune Ki

K_12 = Q11.K * tune1
Ki_12 = Q11.Ki * factor_for_K
Kd_12 = Q11.Kd * tune2
Kp_12 = Q11.Kp * tune3

% K_12 = 1.5
% Kp_12 = 0.23
% Ki_12 = 0.31
% Kd_12 = 5e-5

% Q12.Kp = Kp_12 / factor_for_K;
% Q12.Ki = Ki_12 / factor_for_K;
% Q12.Kd = Kd_12 / factor_for_K;
% Q12.K = K_12 / factor_for_K;

% D_PID = Q12.Kp + Q12.Ki*(1/s) + Q12.Kd*(2*CF*s)/(Q8.Nhat*s+2*CF);
%Q8.Dp = (1/s)*2*CF/(Q8.Nhat*s+2*CF);
%Q12.X = (Q12.K * Q12.D * GH)/ (1 + Q12.K * Q12.D * GH);

% G12 = K_12 * (Kp_12 + Ki_12*1/s + Kd_12*(2*CF*s)/(Q8.Nhat*s+2*CF)) * Q2.Ga * Q4.Gj;
H12 = Q3.Hs * 10^-3 * Q5.Du * Q5.Kfb;
% Q12.X = G12/(1+G12*H12);
% stepinfo(Q12.X,"RiseTimeLimits", [0.1 1])
% 
% figure(); step(Q12.X, 0.05); hold on; grid on


% Q12D = Q12.Kp+Q12.Ki*1/s+Q12.Kd*(2*CF*s)/(Q8.Nhat*s+2*CF);
%%%

Q12.K = K_12*Ki_12;
Q12.Kp = Kp_12/Ki_12;
Q12.Ki = 1;
Q12.Kd = Kd_12/Ki_12;

G12 = Q12.K * (Q12.Kp + Q12.Ki*1/s + Q12.Kd*(2*CF*s)/(Q8.Nhat*s+2*CF)) * Q2.Ga * Q4.Gj;
Q12.X = G12/(1+G12*H12);
step(Q12.X)
stepinfo(Q12.X,"RiseTimeLimits", [0.1 1])
%% Q13
figure(); pzmap(Q12.K * Q12D);
Q13.Z = [-142 -1.31];
[Gm, Q13.PM] = margin(Q12.K*Q12D*Q5.GH);

%% Q14
Bf = 1.5 * Bf;

Q14A= [0 1 0
       -(Kt*Ns^2*Nf^2)/JM -(Br+Bs*Ns^2+(Bf+Bt)*Ns^2*Nf^2)/JM Km/JM
       0 -Km/Lw -Rw/Lw];

phi = inv(s*eye(3)-Q14A);
Q14G = Q4C * phi * Q4B + Q4D;

Q14.Gj = minreal(Q14G(1)*180/pi); 
Q14.G = Q2.Ga * Q14.Gj;

%% Q15
close all
K_15 = Q11.K * 0.8
Ki_15 = Q11.Ki * 0.45;
Kd_15 = Q11.Kd * 0.035;
Kp_15 = Q11.Kp * 1.45;

D_PID = Kp_15 + Ki_15*(1/s) + Kd_15*(2*CF*s)/(Q8.Nhat*s+2*CF);

G15 = K_15 * (Kp_15 + Ki_15*1/s + Kd_15*(2*CF*s)/(Q8.Nhat*s+2*CF)) * Q2.Ga * Q14.Gj;
H15 = Q3.Hs * 10^-3 * Q5.Du * Q5.Kfb;
Q15.X = G15/(1+G15*H15);
stepinfo(Q15.X)
figure(13)
step(Q15.X)

Q15.K = K_15*Ki_15;
Q15.Kp = Kp_15/Ki_15;
Q15.Ki = 1;
Q15.Kd = Kd_15/Ki_15;

Q15.G = Q15.K * (Q15.Kp + Q15.Ki*1/s + Q15.Kd*(2*CF*s)/(Q8.Nhat*s+2*CF)) * Q2.Ga * Q14.Gj;
Q15.X = G15/(1+G15*H15);
step(Q12.X,0.05); hold on; grid on
step(Q15.X,0.05)