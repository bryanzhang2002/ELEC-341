% Solution to ELEC 341 Project 1
% NAME: Bryan Zhang
% SN: 69238335
% DATE: 2022-11-05
% EMAIL: bryan.zhang@alumni.ubc.ca

clear all; close all; clc;

SN = 69238335; sn = [6 9 2 3 8 3 3 5];
A = sn(1)+10; B = sn(2)+10; C = sn(3)+10; D = sn(4)+10; E = sn(5)+10; F = sn(6)+10; G = sn(7)+10; H = sn(8)+10;
s = tf('s');

p1DSPlot(SN);
fig = gcf;
axObjs = fig.Children;
dataObjs = axObjs.Children;
xData = dataObjs(1).XData; 
yData = (dataObjs(1).YData);
clf;
plot(xData*10^(-3), yData, 'k'); grid on; hold on; % plot in seconds
title('Data-Sheet Curve for SN = 69238335'); xlabel('Time(s)'); ylabel('Ouput(V)') 

%% Q1
FV = 13.2;
Tr = 0.000795;
Tp = 0.0016;
Ts = 0.0043;
peak = 16.5;

overshoot = peak - FV;

xline(Tr, '-', 'Tr'); xline(Tp, '-', 'Tp'); xline(Ts, '-', 'Ts');
yline(FV); yline(FV*0.98, '-.'); yline(FV*1.02, '-.', 'FV');

Q1.Tr = Tr*1e3;
Q1.Tp = Tp*1e3;
Q1.Ts = Ts*1e3;
Q1.Pos = overshoot/FV * 100; % Percent Overshoot (%)

%% Q2
zeta = sqrt(log(overshoot/FV)^2/(pi^2+log(overshoot/FV)^2)); beta = sqrt(1-zeta^2);   

% rise/peak time approx
wn_rise = (pi - atan(beta/zeta))/(beta*Q1.Tr/1e3);
wn_peak = pi/(beta*Q1.Tp/1000);

wn = (wn_peak + wn_rise)/2;
Q2.Ga = sysID(FV, zeta, wn); hold on; step(Q2.Ga);

%% Q3
syms U Y X1 s
G1 = 100/(s+A*500);
G2 = 100/(s+B*600);
G3 = 10^5/(s+C*700);
G4 = 1500/(s+D*800);
H1 = 3/(s+F*5);

U = 1;
eqn1 = (U*G1 + X1)*G2 == Y;
eqn2 = (U - (Y+X1*G4)*H1)*G3 == X1;
[X1, Y] = solve(eqn1, eqn2, X1, Y);

[num, den] = numden(Y/U);
Q3G = tf(sym2poly(num), sym2poly(den));

Q3.Hs = minreal(Q3G);
Q3.Kdc = dcgain(Q3G);
figure(); pzmap(Q3.Hs);

clear s; s = tf('s');
%% Q4
% Motor and Lead Screw
Rw = A/2;          % Ω
Lw = B*30   *1e-6; % H
Km = C      *1e-3; % Nm/A

Mh = (D+E)  *1e-3; % kg

Jr = F/15   *1e-7; % kg*m^2
Br = G/30   *1e-6; % Nms

Js = H/5    *1e-7; % kg*m^2
Ms = A/4    *1e-3; % kg
Bs = B/3;          % Ns/m 
Ns = 3e-2/(2*pi);  % m/rad

% Mechanism & Controller
Jf = C/3        *1e-73*3; % kg*m^2
Bf = D/50*3;             % Nms
Nf = 10*pi/180  *1e2;  % rad/m

Bt = E          *1e-3*3; % Nms
Kt = F*30       *1e-3*3; % Nm

L6 = 100        *1e-3; % m

CF = 200;              % Hz
g = 9.81;              % m/s^2

JM = Jr + Js + Ns^2*(Mh+Ms) + Ns^2*Nf^2*Jf;

Q4A= [0 1 0
       -(Kt*Ns^2*Nf^2)/JM -(Br+Bs*Ns^2+(Bf+Bt)*Ns^2*Nf^2)/JM Km/JM
       0 -Km/Lw -Rw/Lw];

Q4B = [0 0 1/Lw]';

% Q4C = [180/pi 0 0
%         Kt*(Ns*Nf)^2 Bt*(Ns*Nf)^2 0];

Q4C = [1 0 0
        Kt*(Ns*Nf)^2/L6 0 0];
Q4D = [0 0]';

phi = inv(s*eye(3)-Q4A);
Q4G = Q4C * phi * Q4B + Q4D;

Q4.Gj = minreal(Q4G(1)*180/pi); 
Q4.Gt = Q4G(2)*1e3/(g*Ns*Nf*3); figure(); step(Q4.Gt);

%% Q5
Q5.Du = CF/(s+CF);
Q5.Kfb = 1/(Q3.Kdc*1e-3);
Q5.GH = (Q2.Ga*Q4.Gj*Q3.Hs*1e-3*Q5.Du*Q5.Kfb);

%% Q6
% Q6.Xj = (Q4.Gj*Q2.Ga/(1+Q5.GH));
G6 = Q2.Ga*Q4.Gj;
H6 = Q3.Hs*1E-3*Q5.Du*Q5.Kfb;
Q6.Xj = G6/(1+G6*H6);

%% Q7
Q7.Xt = minreal(Q4.Gt*Q2.Ga/(1+Q5.GH));

%% Verfication
Q4.Gj_KEY = 1.011e09/(s^3 + 1.901e04* s^2 + 7.006e07* s + 9.559e07);
figure(); step(Q4.Gj); hold on; step(Q4.Gj_KEY); title("Q4.Gj");
Q4.Gt_KEY = 5.847e08/(s^3 + 1.901e04 *s^2 + 7.006e07 *s + 9.559e07);
figure(); step(Q4.Gt); hold on; step(Q4.Gt_KEY); title("Q4.Gt");

Q5.GH_KEY = (1.823e27*s^8 + 1.193e32* s^7 + 3.339e36* s^6 + 5.178e40* s^5 + 4.808e44* s^4 + 2.676e48* s^3 + 8.289e51* s^2 + 1.117e55* s + 6.921e56)/(s^16 + 1.066e05* s^15 + 5.16e09* s^14 + 1.502e14* s^13 + 2.927e18* s^12 + 4.03e22* s^11 + 4.028e26* s^10 + 2.96e30 *s^9 + 1.603e34 *s^8 + 6.353e37* s^7 + 1.82e41* s^6+3.677e44* s^5 + 4.958e47 *s^4 + 3.764e50* s^3 + 7.695e52* s^2 + 3.792e54* s + 5.033e54);
figure(); step(Q5.GH); hold on; step(Q5.GH_KEY); title("Q5.GH")

Q6.Xj_KEY = (8.782e16*s^11 + 7.505e21*s^10 + 2.846e26*s^9 + 6.286e30*s^8 + 8.919e34*s^7 + 8.439e38*s^6 + 5.335e42*s^5 + 2.183e46*s^4 + 5.315e49*s^3 + 6.221e52*s^2+ 1.394e55*s + 6.921e56)/(s^16 + 1.066e05*s^15 + 5.16e09*s^14 + 1.502e14*s^13 + 2.927e18*s^12 + 4.03e22*s^11 + 4.028e26*s^10 + 2.96e30*s^9 + 1.603e34*s^8 + 6.353e37*s^7 + 1.82e41*s^6++ 3.678e44*s^5 + 4.963e47*s^4 + 3.791e50*s^3 + 8.524e52*s^2 + 1.497e55*s + 6.972e56);
figure(); step(Q6.Xj); hold on; step(Q6.Xj_KEY); title("Q6.Xj");

Q7.Xt_KEY = (5.078e16* s^13 + 5.305e21* s^12 + 2.506e26* s^11 + 7.067e30* s^10 + 1.322e35* s^9 + 1.723e39* s^8 + 1.597e43* s^7 + 1.054e47* s^6 + 4.868e50 *s^5 + 1.504e54* s^4+2.844e57* s^3 + 2.673e60* s^2 + 5.722e62* s + 2.803e64)/(s^18 + 1.256e05* s^17 + 7.256e09 *s^16 + 2.557e14* s^15 + 6.143e18* s^14 + 1.065e23* s^13 + 1.374e27* s^12 + 1.344e31* s^11 + 1.005e35* s^10 + 5.755e38* s^9+ 2.512e42* s^8 + 8.276e45* s^7 + 2.023e49 *s^6 + 3.557e52 *s^5 + 4.205e55* s^4 + 2.819e58 *s^3 + 6.255e60 *s^2 + 1.062e63 *s + 4.883e64);
figure(); step(Q7.Xt); hold on; step(Q7.Xt_KEY); title("Q7.Xt");

function [f] = sysID(FV, zeta, wn)
    num = [FV*wn^2];
    den = [1 2*zeta*wn wn^2];
    f = tf(num,den);
end