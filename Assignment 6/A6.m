% Solution to ELEC 341 Assignment 7
% NAME: Bryan Zhang
% SN: 69238335
% DATE: 2022-11-15
% EMAIL: bryan.zhang@alumni.ubc.ca

clear all; close all; clc;

SN = 69238335; sn = [6 9 2 3 8 3 3 5];
A = sn(1)+10; B = sn(2)+10; C = sn(3)+10; D = sn(4)+10; E = sn(5)+10; F = sn(6)+10; G = sn(7)+10; H = sn(8)+10;
s = tf('s');

%% Q1
CF = 10*F; zeroes = -5*C; poles = [-3*B -A (-1+j)*2*D (-1-j)*2*D];
syms k; k_adjusted = sym2poly(solve(k*5*C/(3*B*A*2*(2*D)^2) == 300/G,k));

EMS = zpk(zeroes, poles, k_adjusted);
SEN = 25*E/(s+25*E);
Kf = 1/dcgain(SEN);
Df = CF/(s+CF);

Q1.GH = EMS*SEN*Kf*Df;
Q1.Ku = margin(Q1.GH);

Q1K = Q1.Ku/2;
Q1.X = Q1K*EMS/(1+Q1K*Q1.GH);
Q1.Ess = dcgain(1/(1+Q1K*Q1.GH))*1e2;

%% Q2
Q2.GH = 1/s*EMS*SEN*Kf*Df;
Q2.Ku = margin(Q2.GH);
Q2K = Q2.Ku/2;
Q2.X = Q2K*EMS*(1/s)/(1+Q2K*(1/s)*Q1.GH);
Q2.Ess = dcgain(1/(1+Q2K*(1/s)*Q1.GH));

%% COW
% figure(); pzmap(Q1.GH); title("pole-zero map of open-loop transfer function for velocity");
% figure(); pzmap(Q2.GH); title("pole-zero map of open-loop transfer function for position");
% figure(); impulse(Q1.GH); hold on; impulse(Q2.GH); impulse(Q1.X); impulse(Q2.X); title("impulse responses of velocity and position OL/CL TFs"); legend('velocity OL', 'displacement OL', 'velocity CL', 'displacement CL'); grid on; hold off;
% figure(); step(Q1.GH); hold on; step(Q2.GH); step(Q1.X); step(Q2.X); title("step responses of velocity and position OL/CL TFs"); legend('velocity OL', 'displacement OL', 'velocity CL', 'displacement CL'); grid on; hold off;

%% Q3
%dcgain(1/(1+Q3.K1 * Q1.GH)) = 0.3
Q3.K1 = (1/0.3 - 1)/dcgain(Q1.GH);

index = 1;
k_vec = [Q1K:-0.0001:0];
overshoot_vec = zeros(1, length(k_vec));

for k3 = k_vec
    overshoot_vec(index) = stepinfo(k3*EMS/(1+k3*Q1.GH)).Overshoot;
    index = index + 1;
end

Q3K = k_vec(overshoot_vec == 0);
Q3.K2 = Q3K(1); % Q3.K2 = 0.0081;
step(Q1.X); hold on; step(Q3.K2*EMS/(1+Q3.K2*Q1.GH)); grid on;

%% Q4
% 10 percent overshoot
index = 1;
k_vec = [Q2K:-0.0001:0];
overshoot_vec = zeros(1, length(k_vec));

for k4 = k_vec
    overshoot_vec(index) = stepinfo(k4*EMS*(1/s)/(1+k4*(1/s)*Q1.GH)).Peak - 1;
    index = index + 1;
end

Q4K = k_vec(overshoot_vec<0.1);
Q4.K1 = Q4K(1);% Q4.K1 = 0.2376;
figure(); step(Q2.X); hold on; step(Q4.K1*EMS*(1/s)/(1+Q4.K1*(1/s)*Q1.GH));

% 0 percent overshoot
index = 1;
k_vec = [Q2K:-0.0001:0];
overshoot_vec = zeros(1, length(k_vec));

for k4 = k_vec
    overshoot_vec(index) = stepinfo(k4*EMS*(1/s)/(1+k4*(1/s)*Q1.GH)).Peak - 1;
    index = index + 1;
end

Q4K = k_vec(overshoot_vec < 1e-6);
Q4.K2 = Q4K(1); % Q4.K2 = 0.1260;
step(Q4.K2*EMS*(1/s)/(1+Q4.K2*(1/s)*Q1.GH)); grid on; yline(1.1);