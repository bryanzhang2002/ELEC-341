%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution to Assignment 8
% 2022/11/18
% Katie Seifert
% Student #68469311
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clean up workspace
clear all; clc;

% SN variable must contain Student Number
% This must be right or solution will not be graded
SN    = 68469311;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Variables

A = 6 + 10;
B = 9 + 10;
C = 2 + 10;
D = 3 + 10;
E = 8 + 10;
F = 3 + 10;
G = 3 + 10;
H = 5 + 10;

s = tf('s');
%Variables
CF = 10*F;
G_DCG = 300/G;
PA = -A;
PB = -3*B;
ZC = -5*C;
PD1 = 2*D*(-1+i);
PD2 = 2*D*(-1-i);
PE = -25*E;

%new derivates
P = -2*CF;


%% Asn6
%Q1
% create a transfer function from the figure 2

EMS = zpk(ZC, [PA, PB, PD1, PD2], 1);
SEN = zpk([], PE, 1);
kf = 1/dcgain(SEN);
k_EMS = 1/dcgain(EMS);
dh = CF/(s+CF);

Q1G = EMS*300/G*k_EMS;
Q1GH = EMS*300/G*SEN*dh*kf*k_EMS;
%pzmap(Q1GH)
%Q1.Ku = margin(Q1.GH);
%Q1.K = Q1.Ku/2;
%Q1.X = Q1.K*Q1.G/(1+Q1.K*Q1.GH); 
%Q1.Ess = dcgain(1/(1+Q1.K*Q1.GH))*100;

%Q2
% same transfer functions

Q2G = Q1G/s;
Q2GH = Q1GH/s;
%pzmap(Q2GH)
%Q2Ku = margin(Q2GH);
%Q2K = Q2Ku/2;
%Q2X = Q2K*Q2G/(1+Q2K*Q2GH);
%Q2Ess = dcgain(1/(1+Q2K*Q2GH))*100;

%% Q1
%pzmap(Q1GH)
Z = PA;
Kp = 1/(-Z)
Dd = (s-Z)/(s*(-Z))
Q1.Kp = Kp;
Q1.D = Dd;

%% Q2
Ku = margin(Dd*Q1GH)
t = 0;
pt = -10*16;
min = 1000;
kf = 0
zf = 0
while t < Ku
    Kt = t;
    while pt < -1
    Z = pt;
    Kp = 1/(-Z);
    Dd = (s-Z)/(s*(-Z));
    X = Dd*Kt*Q1G/(1+Dd*Kt*Q1GH);
    S = stepinfo(X);
    if ( S.Overshoot <= 10)
        disp("/n")
        disp(S.Overshoot)
        disp(Kt)
        disp(Z)
        disp(S.SettlingTime)
        if (min > S.SettlingTime)
            min = S.SettlingTime;
            kf = Kt;
            zf = Z;
            disp("min=")
            disp(min)
        end
    end
    pt = pt + 0.1*16;
    end
    pt = -10*16;
    t = t+0.01*Ku;
end
disp("final values")
disp(min)
disp(kf)
disp(zf)
%% Q2 Final answers
Kp = 1/(-zf);
Dd = (s-zf)/(s*(-zf));
X = Dd*kf*Q1G/(1+Dd*kf*Q1GH);
Q2.K = kf;
Q2.Z = zf;
Q2.X = minreal(Dd*kf*Q1G/(1+Dd*kf*Q1GH));

%% Q3
%D = -p/(z1*z2)*(s-z1)*(s-z2)/s(s-p)
%Kd = (Z-P)/(P*Z)
%Kp = 1/(-Z)

Z1 = PA
Z2 = PB

Dd = -P/(Z1*Z2)*(s-Z1)*(s-Z2)/(s*(s-P))
Kp = 1/P-(Z1+Z2)/(Z1*Z2)
Kd = 1/(Z1*Z2)+Kp/P
Q3.Kp = Kp;
Q3.Kd = Kd;
Q3.D = Dd;
%% Q4
Ku = margin(Dd*Q1GH)
K = 0.5*Ku
X = minreal(Dd*K*Q1G/(1+Dd*K*Q1GH))
Q4.Ku = Ku;
Q4.X = X;

%% Q5
Z1 = PD1
Z2 = PD2

Dd = -P/(Z1*Z2)*(s-Z1)*(s-Z2)/(s*(s-P));
Kp = 1/P-(Z1+Z2)/(Z1*Z2);
Kd = 1/(Z1*Z2)+Kp/P;
Ku = margin(Dd*Q1GH)
K = 0.5*Ku
X = minreal(Dd*K*Q1G/(1+Dd*K*Q1GH))
Q5.Ku = Ku;
Q5.X = X;

%% Q6
%{
min = 1000;
kf = 0;
z1f = 0;
z2f = 0;
for Kt=0:0.01*Ku:Ku
    for zt1 = 10*PA:-0.1*PA:-1
        for zt2 = 10*PB:-0.1*PB:-1
            Kp = 1/(-zt1);
            Kd = P/zt2;
            Dd = -P/(zt1*zt2)*(s-zt1)*(s-zt2)/(s*(s-P));
            X = Dd*Kt*Q1G/(1+Dd*Kt*Q1GH);
            S = stepinfo(X);
            if ( S.Overshoot <= 10)
                disp("/n")
                disp(S.Overshoot)
                disp(Kt)
                disp(Z)
                disp(S.SettlingTime)
                if (min > S.SettlingTime)
                    min = S.SettlingTime;
                    kf = Kt;
                    z1f = zt1;
                    z2f = zt2;
                    disp("min=")
                    disp(min)
                end
            end
        end
    end
end
%}
min = 1000;
kf = 0;
zf = 0;
for Kt=0:0.01*Ku:Ku
    Ztr = PA
    Zti = 0
    while Ztr+Zti > PD1
    %for Zt = 10*PA:-0.1*PA:-1
        Zt1 = Ztr+Zti
        Zt2 = Ztr-Zti
        Dd = -P/(Zt1*Zt2)*(s-Zt1)*(s-Zt2)/(s*(s-P));
        Kp = 1/P-(Zt1+Zt2)/(Zt1*Zt2);
        Kd = 1/(Zt1*Zt2)+Kp/P;
        X = Dd*Kt*Q1G/(1+Dd*Kt*Q1GH);
        S = stepinfo(X);
        if ( S.Overshoot <= 10)
            disp("/n")
            disp(S.Overshoot)
            disp(Kt)
            disp(Zt1)
            disp(Zt2)
            disp(S.SettlingTime)
            if (min > S.SettlingTime)
                min = S.SettlingTime;
                kf = Kt;
                zf1 = Zt1;
                zf2 = Zt2;
                disp("min=")
                disp(min)
            end
        end
        Ztr = Ztr + 0.1*PA;
        Zti = Zti + 0.1*i;
    end
end
disp("final values")
disp(min)
disp(kf)
disp(zf1)
disp(zf2)

%% Q6 Final answers
Dd = -P/(zf1*zf2)*(s-zf1)*(s-zf2)/(s*(s-P))
Kp = 1/P-(zf1+zf2)/(zf1*zf2)
Kd = 1/(zf1*zf2)+Kp/P
X = Dd*kf*Q1G/(1+Dd*kf*Q1GH)
zf1
zf2
Q6.K = kf;
Q6.Z = [zf1 zf2];
Q6.X = minreal(X);