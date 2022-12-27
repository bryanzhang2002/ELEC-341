% Solution to ELEC 341 Assignment 5
% NAME: Bryan Zhang
% SN: 69238335
% DATE: 2022-10-04
% EMAIL: bryan.zhang@alumni.ubc.ca

clear all; close all; clc;

SN = 69238335; sn = [6 9 2 3 8 3 3 5];
A = sn(1)+10; B = sn(2)+10; C = sn(3)+10; D = sn(4)+10; E = sn(5)+10; F = sn(6)+10; G = sn(7)+10; H = sn(8)+10;
s = tf('s');

% Mechanical Circuit
M0 = A; M3 = A; M1 = B; M2 = B; B0 = C; B20 = E; B31 = F; K1 = G; K32 = G; K20 = H; K21 = H; F0 = 300; B1 = D;

% Linear Actuator
Rw = A/3; Lw = B*1e-3; Km = G/2; Vin = 150; M03 = A+B; B03 = C+D; N0=1/2; g = 9.81;

%% Q1
Q1.A = [0 0 0 0 1 0 0 0
        0 0 0 0 0 1 0 0
        0 0 0 0 0 0 1 0
        0 0 0 0 0 0 0 1
        -K20/M0 0 K20/M0 0 -(B0+B20)/M0 0 B20/M0 0
        0 -(K1+K21)/M1 K21/M1 0 0 -B31/M1 0 B31/M1
        K20/M2 K21/M2 -(K20+K21+K32)/M2 K32/M2 B20/M2 0 -B20/M2 0
        0 0 K32/M3 -K32/M3 0 B31/M3 0 -B31/M3];

Q1.B = [0 0 0 0 F0/M0 0 0 0]';

%% Q2
Q2.C = [0 0 0 1 0 0 0 0; K20 0 -K20 0 0 0 0 0];
Q2.D = [0;0];

% phi = inv(s*eye(8)-Q1.A);
% Q1.G = Q2.C * phi * Q1.B + Q2.D;
% impulse(Q1.G(1));

%% Q3
Q3A =  [0 0 0 0 1 0 0 0 0
        0 0 0 0 0 1 0 0 0
        0 0 0 0 0 0 1 0 0
        0 0 0 0 0 0 0 1 0
        -K20/M03 0 K20/M03 0 -(B03+B20)/M03 0 B20/M03 0 Km/(2*M03)
        0 -(K1+K21)/M1 K21/M1 0 0 -B31/M1 0 B31/M1 0
        K20/M2 K21/M2 -(K20+K21+K32)/M2 K32/M2 B20/M2 0 -B20/M2 0 0
        0 0 K32/M3 -K32/M3 0 B31/M3 0 -B31/M3 0
        0 0 0 0 -Km/Lw 0 0 0 -Rw/Lw];

Q3B = [0 0 0 0 0 0 0 0 1/Lw]';
Q3C = [0 0 0 1 0 0 0 0 0];
Q3D = 0;

phi3 = inv(s*eye(9)-Q3A);
Q3.G = Q3C * phi3 * Q3B + Q3D;
% impulse(Q1.G(1));

%% Q4
Q4B = [0 0 0 0 0 0 0 0 1/Lw; 0 0 0 0 -g -g -g -g 0]';
Q4D = [0 0];

phi4 = inv(s*eye(9) - Q3A);
Q4G = Q3C * phi4 * Q4B + Q4D;

Q4.Gv = Q4G(1);
Q4.Gg = Q4G(2);