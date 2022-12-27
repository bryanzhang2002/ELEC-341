% Solution to ELEC 341 Assignment 4
% NAME: Bryan Zhang
% SN: 69238335
% DATE: 2022-10-04
% EMAIL: bryan.zhang@alumni.ubc.ca

clear all; clc; hold on;

SN = 69238335; sn = [6 9 2 3 8 3 3 5];
A = sn(1)+10; B = sn(2)+10; C = sn(3)+10; D = sn(4)+10; E = sn(5)+10; F = sn(6)+10; G = sn(7)+10; H = sn(8)+10;
syms s U1 U2 Y1 Y2

% Block Diagram model
G1 = 1/(s+A); G2 = 10/(s+B); G3 = 10/(s+C); G4 = 10/(s+D); G5 = 1/(s+E); H1 = 100/(s+F); H2 = 1000/(s+G);

% Motor Model
Rw = 1 + A/10; Lw = (100 + 10*B)*1e-6; Jr = (C/10)*1e-6; Br = (D+E+F)*1e-6; Km = (10+G)*1e-3; Vin = 12; Jf = (G/30)*1e-6; Bf = H*Br;

% Q1
U2 = 0;
eqn1 = (U1 - Y2*H2)*G1*G2 == Y1;
eqn2 = ((U1 - Y2*H2)*G1*G3 + (Y1*H1-U2)*G4)*G5 == Y2;

[Y1, Y2] = solve(eqn1, eqn2, Y1, Y2);
[num11, den11] = numden(Y1/U1);
Q1.G11 = minreal(tf(sym2poly(num11), sym2poly(den11)));
% [y11, t11] = step(Q1.G11); plot(t11, y11, 'r', 'LineWidth', 2);

[num12, den12] = numden(Y2/U1);
Q1.G12 = minreal(tf(sym2poly(num12), sym2poly(den12)));
% [y12, t12] = step(Q1.G12); plot(t12, y12, 'r', 'LineWidth', 2);

clear U2 Y1 Y2; syms U2 Y1 Y2; U1 = 0;
eqn1 = (U1 - Y2*H2)*G1*G2 == Y1;
eqn2 = ((U1 - Y2*H2)*G1*G3 + (Y1*H1-U2)*G4)*G5 == Y2;

[Y1, Y2] = solve(eqn1, eqn2, Y1, Y2);
[num21, den21] = numden(Y1/U2);
Q1.G21 = minreal(tf(sym2poly(num21), sym2poly(den21)));
% [y21, t21] = step(Q1.G21); plot(t21, y21, 'r', 'LineWidth', 2);

[num22, den22] = numden(Y2/U2);
Q1.G22 = minreal(tf(sym2poly(num22), sym2poly(den22)));
% [y22, t22] = step(Q1.G22); plot(t22, y22, 'r', 'LineWidth', 2);

% Q2
clear U2 Y1 Y2; syms U Y1 Y2;
eqn1 = (U - Y2*H2)*G1*G2 == Y1;
eqn2 = ((U - Y2*H2)*G1*G3 + (Y1*H1-U)*G4)*G5 == Y2;

[Y1, Y2] = solve(eqn1, eqn2, Y1, Y2);
[num1, den1] = numden(Y1/U);
Q2.G1 = minreal(tf(sym2poly(num1), sym2poly(den1)));
% [y1, t1] = step(Q2.G1); plot(t1, y1, 'r', 'LineWidth', 2);

[num2, den2] = numden(Y2/U);
Q2.G2 = minreal(tf(sym2poly(num2), sym2poly(den2)));
% [y2, t2] = step(Q2.G2); plot(t2, y2, 'r', 'LineWidth', 2);

Q2.G1 = Q1.G11 + Q1.G21;
% step(Q2.G1)
Q2.G2 = Q1.G12 + Q1.G22;
% step(Q2.G2)

% Q3
s = tf('s');
Jtot = Jr + Jf;
Btot = Br + Bf;

Q3.Ye = minreal(1/(s*Lw + Rw));
Q3.Ym = minreal(1/(s*Jtot+Btot));
Q3.G = minreal(feedback(Q3.Ye*Km*Q3.Ym, Km));