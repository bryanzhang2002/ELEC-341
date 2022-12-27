% Solution to ELEC 341 Assignment 1
% NAME: Bryan Zhang
% SN: 69238335
% DATE: 2022-09-17
% EMAIL: bryan.zhang@alumni.ubc.ca

clear all; close all; clc;

SN = 69238335; sn = [6 9 2 3 8 3 3 5];
A = sn(1)+10; B = sn(2)+10; C = sn(3)+10; D = sn(4)+10; E = sn(5)+10; F = sn(6)+10; G = sn(7)+10; H = sn(8)+10;
s = tf('s');

% Circuit Parameters
R1 = 10 * A; L1 = B * 1e-3; C1 = C * 1e-6; R2 = 10 * D; L2 = E * 1e-3; C2 = F * 1e-6; Vin = G; Iin = H;

%% Q1
syms s Vo1 Vout
KCL1 = Vin/RR(1/(s*C1), R1) == -Vo1/R1 + -Vo1 * s*C1 - Vo1/(s*L1);
KCL2 =  Vo1/(s*L2) + Vo1/R2 == -Vout/RR(1/(s*C2), R2);

[Vo1, Vout] = solve(KCL1, KCL2, Vo1, Vout);
[num1, den1] = numden(Vout/Vin);
Q1.G = minreal(tf(sym2poly(num1), sym2poly(den1)));

%% Q2
step(Vin * Q1.G); grid on; hold on;
% from https://www.dcode.fr/function-equation-finder
x = 0:0.001:0.06; y1 = 11.1953 + 13.8364 * exp(-278.424*x);
plot(x,y1, 'r--');

Q2.K = 13.8364;
Q2.A = 278.424;

%% Q3
figure();
[y, t] = impulse(Q1.G);
plot (t, y);

Q3.Tl = (t(2) - t(1));
Q3.Vh = 1/Q3.Tl;
u = Q3.Vh*(heaviside(t)-heaviside(t-Q3.Tl));

[ys, ts] = lsim(Q1.G, u, t);
hold on; plot(ts, ys); grid on;
Q3.Vh = Q3.Vh * 1e-3;

%% Q4
syms v1 v2 v3
KCL1 = Iin == (v1-v3)/(R2+s*L2+1/(s*C2)) + (v1-v2)/R2;
KCL2 = (v1-v2)/R2 == (v2-v3)/(s*L2) + v2/RR(s*L1, R1+1/(s*C1));
KCL3 = (v1-v3)/(R2+s*L2+1/(s*C2))+(v2-v3)/(s*L2) == v3*s*C2;

[v1, v2, v3] = solve(KCL1, KCL2, KCL3, v1, v2, v3);
Iout = v2/RR(s*L1, R1+1/(s*C1));

[num4, den4] = numden(Iout/Iin);
Q4.G = minreal(tf(sym2poly(num4), sym2poly(den4)));

%% Q5
figure(); step(Iin*Q4.G); grid on;
% from https://www.dcode.fr/function-equation-finder
x = 0:0.001:0.06; y1 = 15.0054 + 11.3839 * exp(-90.9168*x);
hold on; plot(x,y1,'r--');
Q5.K = 11.3839;
Q5.A = 90.9168;

%% Q6
tau = 1/Q5.A;
Q6.Ts = 5 * tau * 1e3;


function [Z] = RR(Z1, Z2)
    Z = (Z1 * Z2)/(Z1+Z2);
end