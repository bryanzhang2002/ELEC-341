% Solution to ELEC 341 Assignment 7
% NAME: Bryan Zhang
% SN: 69238335
% DATE: 2022-11-04
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

% Q1K = Q1.Ku/2;
% Q1.X = Q1K*EMS/(1+Q1K*Q1.GH);
% Q1.Ess = dcgain(1/(1+Q1K*Q1.GH))*1e2;

% add the lead controller dynamics to this (Add a pole at 2CF, use a zero
% to cancel out right most pole)
p = 2*CF;
z = A;

% Q1.Kd = p-z/(p*z);
% Q1.D = 1+(p-z)/(z*p)*(s*p)/(s+p);
Q1.Kd = (p-z)/(p*z); % CORRECT
Q1.D = 1+(p-z)/(z*p)*(s*p)/(s+p); pzmap(Q1.D*Q1.GH)

%% Q2
Ku = margin(Q1.D*Q1.GH);
K = Ku/2;
Q2X = (K*Q1.D*EMS)/(1+K*Q1.D*Q1.GH);
Q2.Ess = dcgain(1/(1+K*Q1.D*Q1.GH))*1e2;
step(K*Q1.D*EMS/(1+K*Q1.D*Q1.GH))

K = 0.99*Ku;
Q2X = (K*Q1.D*EMS)/(1+K*Q1.D*Q1.GH);
Q2.Ess99 = dcgain(1/(1+K*Q1.D*Q1.GH))*1e2;

%% Q2
% Q2.GH = 1/s*EMS*SEN*Kf*Df;
% Q2.Ku = margin(Q2.GH);
% Q2K = Q2.Ku/2;
% Q2.X = Q2K*EMS*(1/s)/(1+Q2K*(1/s)*Q1.GH);
% Q2.Ess = dcgain(1/(1+Q2K*(1/s)*Q1.GH));

%% COW
% figure(); pzmap(Q1.GH); title("pole-zero map of open-loop transfer function for velocity");
% figure(); pzmap(Q2.GH); title("pole-zero map of open-loop transfer function for position");
% figure(); impulse(Q1.GH); hold on; impulse(Q2.GH); impulse(Q1.X); impulse(Q2.X); title("impulse responses of velocity and position OL/CL TFs"); legend('velocity OL', 'displacement OL', 'velocity CL', 'displacement CL'); grid on; hold off;
% figure(); step(Q1.GH); hold on; step(Q2.GH); step(Q1.X); step(Q2.X); title("step responses of velocity and position OL/CL TFs"); legend('velocity OL', 'displacement OL', 'velocity CL', 'displacement CL'); grid on; hold off;

%% Q3
% Requirments: Peak value < 1.2
%              Ess < Q2.Ess
% Goals: Ess as small as possible

KValues = 0.01:0.01:Ku; % Start at 0.01 because 0 gives weird values
ZValues = 1:1:p; % Start at 1 because 0 gives weird values
PValues = zeros(length(ZValues), length(KValues));
EssValues = zeros(length(ZValues), length(KValues));

for ii = 1:length(ZValues)
    Q3D =(s*ZValues(ii))/(s+2*CF);
    Ku=margin(Q3D*GH);
    for jj = 1:length(KValues)
        X = KValues(jj)*D*EMStf/(1+KValues(jj)*GH*D);
        PValues(ii, jj) = stepinfo(X).Peak;
        EssValues(ii, jj) = 1/(1+dcgain(KValues(jj)*D*GH*D))*100;
    end
    fprintf("ii= %d\n", ii)
end
acceptableP = PValues < 1.2;
invert = ~acceptableP;
EssValues(invert) = inf;
[x,y] = find(EssValues == min(min(EssValues)));
Q3.K = KValues(y);
Q3.Z = ZValues(x);
D = (s+Q3.K)/(s+2*CF);
Q3.X = Q3.K*D*EMS/(1+Q3.K*Q1.GH);


%% Q4
Q4.Ku = margin(GH*Q1.D/s);
Q4.X = Q4.Ku/2*Q1.D*EMStf/(s*(1+Q4.Ku/2*Q1.D*GH));