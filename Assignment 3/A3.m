% Solution to ELEC 341 Assignment 3
% NAME: Bryan Zhang
% SN: 69238335
% DATE: 2022-10-04
% EMAIL: bryan.zhang@alumni.ubc.ca

% clear all; close all; clc;

SN = 69238335; sn = [6 9 2 3 8 3 3 5];
A = sn(1)+10; B = sn(2)+10; C = sn(3)+10; D = sn(4)+10; E = sn(5)+10; F = sn(6)+10; G = sn(7)+10; H = sn(8)+10;
s = tf('s');

% Circuit Parameters
M0 = A; M3 = A; M1 = B; M2 = B; B0 = C; B20 = E; B31 = F; K1 = G; K32 = G; K20 = H; K21 = H; F0 = 300; B1 = D;
% Mechanical Circuit
R0 = 1/B0; R20 = 1/B20; R31 = 1/B31; L20 = 1/K20; L32 = 1/K32; L1 = 1/K1; L21 = 1/K21; I0 = F0/s; C0 = M0; C1 = M1; C2 = M2; C3 = M3; R1 = 1/B1;
% Motor Model
Rw = 1 + A/10; Lw = (100 + 10*B)*1e-6; Jr = (C/10)*1e-6; Br = (D+E+F)*1e-6; Km = (10+G)*1e-3; Vin = 12;  

%% Q1
Y = [1/R0+1/R20+s*C0+1/(s*L20) 0 -1/R20-1/(s*L20) 0
    0 1/(s*L1)+s*C1+1/(s*L21)+1/R31 -1/(s*L21) -1/R31
    -1/R20-1/(s*L20) -1/(s*L21) 1/(s*L21)+1/(s*L32)+s*C2+1/(s*L20)+1/R20 -1/(s*L32)
    0 -1/R31 -1/(s*L32) 1/R31+s*C3+1/(s*L32)];
I = [F0; 0; 0; 0];
V = Y\I;

Q1.G = minreal((1/s)*(V(4)/F0)); impulse(Q1.G); grid on;

%% Q2
Q2.G = minreal((V(1) - V(3))/(s*L20*I0)); figure(); impulse(Q2.G); grid on;

%% Q3
Y3 = [1/R0+1/R20+s*C0+1/(s*L20) 0 -1/R20-1/(s*L20) 0
    0 s*C1+1/R1+1/(s*L21)+1/R31 -1/(s*L21) -1/R31
    -1/R20-1/(s*L20) -1/(s*L21) 1/(s*L21)+1/(s*L32)+s*C2+1/(s*L20)+1/R20 -1/(s*L32)
    0 -1/R31 -1/(s*L32) 1/R31+s*C3+1/(s*L32)];
I3 = [F0; 0; 0; 0];
V3 = Y3\I3;

Q3.G = minreal((1/s)*(V3(4)/(F0))); figure(); step(Q3.G); grid on;

%% Q4
Q4.G = minreal(Km/((Lw*Jr)*s^2+(Lw*Br+Jr*Rw)*s+Rw*Br+Km^2));

%% Q5
Q5.G = minreal((Br+s*Jr)/((Lw*Jr)*s^2+(Lw*Br+Jr*Rw)*s+Rw*Br+Km^2));