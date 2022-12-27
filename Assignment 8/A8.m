% Solution to ELEC 341 Assignment 8
% NAME: Bryan Zhang
% SN: 69238335
% DATE: 2022-11-21
% EMAIL: bryan.zhang@alumni.ubc.ca

close all; clear; clc;

SN = 69238335; sn = [6 9 2 3 8 3 3 5];
A = sn(1)+10; B = sn(2)+10; C = sn(3)+10; D = sn(4)+10; E = sn(5)+10; F = sn(6)+10; G = sn(7)+10; H = sn(8)+10;
s = tf('s');


CF = 10*F; zeroes = -5*C; poles = [-3*B -A (-1+1j)*2*D (-1-1j)*2*D];
syms k; k_adjusted = sym2poly(solve(k*5*C/(3*B*A*2*(2*D)^2) == 300/G,k));   % find the k that ensures a gain of 300/G

EMS = zpk(zeroes, poles, k_adjusted);
SEN = 25*E/(s+25*E);
Kf = 1/dcgain(SEN);
Df = CF/(s+CF);

Q1.GH = EMS*SEN*Kf*Df;

%% Q1
% add the LAG controller dynamics to this (Add a pole @ 0 and zero @ most
% dominant pole)
z = A;
Q1.D = 1/z*(s+z)/s;
Q1.Kp = 1/z;

Ku = margin(Q1.D*Q1.GH);
K = Ku/2;
Q1.X = K*Q1.D*EMS/(1+K*Q1.D*Q1.GH);
%COW
figure(); pzmap(Q1.D*Q1.GH); title("Q1: pole @ 0, zero @ A")

%% Q2
% Requirments:
%   overshoot <= 10%
%   Ts as small as possible   

t = 0;          % initial value of master gain
pt = -10*A;     % initial position of zero (negaive convention)
min = 1000;     % initial settling time
kf = 0;         % stores the final value of master gain
zf = 0;         % stores the final value of zero position

%tune the master gain, and then tune all possible values of the zero from -10*A to -1. 
while t < Ku
    Kt = t;
    while pt < -1   % USING NEGATIVE CONVENTION HERE
        Z = pt; 
        Dd = (s-Z)/(s*(-Z));    % USING NEGATIVE CONVENTION HERE
        X = Dd*Kt*EMS/(1+Dd*Kt*Q1.GH);

        if (stepinfo(X).Overshoot <= 10)
            disp(stepinfo(X).Overshoot)                         
            disp(Kt)
            disp(Z)
            disp(stepinfo(X).SettlingTime)
            if (min > stepinfo(X).SettlingTime)
                min = stepinfo(X).SettlingTime;
                kf = Kt;
                zf = Z;
                disp("min=")
                disp(min)
            end
        end
        pt = pt + 0.1*A;
    end
    pt = -10*A;
    t = t+0.01*Ku;
end
disp("final values")
disp(min)
disp(kf)
disp(zf)

Q2.K = kf;
Q2.Z = zf;
Q2.D = 1/-Q2.Z*(s-Q2.Z)/s;      % NEGATIVE CONVENTION
Q2.X = minreal(Q2.D*Q2.K*EMS/(1+Q2.D*Q2.K*Q1.GH)); 

% COW
step(Q1.X); grid on; hold on; step(Q2.X); title("Q1.X vs Q2.X")

%% Q3
% add the PID controller dynamics to this (pole @ [-2CF, 0] zero
% @ two most dominant poles)
p = -2*CF; z1 = -A; z2 = -3*B; Ki = 1;

Q3.Kp = 1/p - (z1+z2)/(z1*z2);
Q3.Kd = 1/(z1*z2) + Q3.Kp/p;
Q3.D = Q3.Kp + Ki*(1/s) + Q3.Kd*(-p*s)/(s - p);
figure(); pzmap(Q3.D*Q1.GH); title("Q3: pole @ [-2CF 0], zero @ two most dominant poles"); grid on;

%% Q4
Q4.Ku = margin(Q3.D*Q1.GH);
K4 = Q4.Ku / 2;
Q4.X = K4*Q3.D*EMS/(1+K4*Q3.D*Q1.GH); 

%% Q5
z1 = (-1+1j)*2*D; z2 = (-1-1j)*2*D;
Q5.Kp = 1/p - (z1+z2)/(z1*z2);
Q5.Kd = 1/(z1*z2) + Q5.Kp/p;
Q5.D = (Q5.Kp + Ki*(1/s)) + Q5.Kd*(-p*s)/(s-p);

Q5.Ku = margin(Q5.D*Q1.GH); K5 = Q5.Ku/2;
Q5.X = K5*Q5.D*EMS/(1+K5*Q5.D*Q1.GH);

% COW
figure(); pzmap(Q5.D*Q1.GH); title("Q5: pole @ [-2CF, 0], zero @ complex poles"); grid on;
figure(); step(Q4.X); hold on; step(Q5.X); grid on; title ("Q4.X vs Q5.X")

%% Q6
% Requirments:
% overshoot <= 10%
% Ts as small as possible
% P = -2*CF; 
min = 1000;
kf = 0;
zf = 0;
for Kt=0:0.01*Q5.Ku:Q5.Ku
    Ztr = -A;
    Zti = 0;
    while Ztr+Zti > 2*D*(-1+1j)
        Zt1 = Ztr+Zti;
        Zt2 = Ztr-Zti;
        Dd = -p/(Zt1*Zt2)*(s-Zt1)*(s-Zt2)/(s*(s-p));
        Kp = 1/p-(Zt1+Zt2)/(Zt1*Zt2);
        Kd = 1/(Zt1*Zt2)+Kp/p;
        X = Dd*Kt*EMS/(1+Dd*Kt*Q1.GH);
        if (stepinfo(X).Overshoot <= 10)
            disp("/n")
            disp(stepinfo(X).Overshoot)
            disp(Kt)
            disp(Zt1)
            disp(Zt2)
            disp(stepinfo(X).SettlingTime)
            if (min > stepinfo(X).SettlingTime)
                min = stepinfo(X).SettlingTime;
                kf = Kt;
                zf1 = Zt1;
                zf2 = Zt2;
                disp("min=")
                disp(min)
            end
        end
        Ztr = Ztr + 0.1*-A;
        Zti = Zti + 0.1*1j;
    end
end
disp("final values")
disp(min)
disp(kf)
disp(zf1)
disp(zf2)

Q6.K = kf;
Q6.D = -p/(zf1*zf2)*(s-zf1)*(s-zf2)/(s*(s-p));
Q6.Z = [zf1 zf2];
Q6.X = minreal(Q6.K*Q6.D*EMS/(1+Q6.K*Q6.D*Q1.GH));

% COW
figure(); step(Q1.X); hold on; step(Q6.X); grid on; title ("Q1.X vs Q6.X")