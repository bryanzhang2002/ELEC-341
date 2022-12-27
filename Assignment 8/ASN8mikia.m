%ASN 8
%Mikia Whitehead
%54514740

clear; clc;

SN = 54514740;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = tf('s');

A = 15; 
B = 14; 
C = 15; 
D = 11; 
E = 14;
F = 17; 
G = 14;
H = 10;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CF = 10*F;
SEN_gain = 1;
EMS_gain = 300/G;

Df = CF/(s + CF);

EMS_zeros = (s + (5*C));
EMS_poles = (s + A) * (s+(3*B)) * (s+(2*D)+(2*D)*1i)*(s+(2*D)-(2*D)*1i);
SEN_poles = (s+(25*E));

SEN = (25*E)/SEN_poles;
EMS_xfer = EMS_zeros/EMS_poles;

Kf = 1/dcgain(SEN);
w = dcgain(EMS_xfer);

EMS = (1/(w/(EMS_gain)))*EMS_xfer;

GH = EMS* SEN * Kf * Df;
%Q1.Ku = margin(open_loop_xfer);
%K = Q1.Ku/2;

%Q1.Ess = dcgain(1/(1 + K*open_loop_xfer))*100;


%% Question 1

% compute derivitage gain Kp
% compute the controller dynamics D
% Compute the ultimate gain of DGH (Ku)
% loop up the settle time of the closed-loop xfer function when K = Ku*50%

%replace the lead controllerr (KD) from the previous assignment 
%with the lag controller (KD)

% controller at pole 0
%use controller zero to cancel the most dominant real system pole

z = -A;
p = 0;
Ki = 1;

Q1.Kp = -1/z;;

%NumD = [Q1.Kp-p*Kd -Kp*p+Ki -Ki*p];
%DenD = [1 -p];

Q1.D = Q1.Kp + Ki * (1/s);

DGH = GH * Q1.D;

Q1.Ku = margin(DGH);

%% Question 2

%Kp moves to zero
%K moves the poles along the root locus

%overshoot <= 10%
% Ts as small as possible

%Ku = margin(open_loop_xfer*Q1.D); 

Q2.K = Q1.Ku * 0.2; %%%%%%%%%%%%%%%%%%%%%%%%%%adjust this for design
Q2.Z = -A; %ajust this for design

D2 = (s-Q2.Z)/(s*-Q2.Z);
DGH2 = D2*GH;
Q2.X = (Q2.K * EMS * D2)/(1 + Q2.K * DGH2);
stepinfo(Q2.X)

%% Question 3

p3 = -2*CF;
z1 = -A;
z2 = -3*B;

Q3.Kp = (1/p3) - (z1+z2)/(z1*z2);
Q3.Kd = 1/(z1*z2) + (Q3.Kp/p3);
Q3.D = Q3.Kp + Ki*(1/s) + Q3.Kd*(-p3*s)/(s - p3);


%% Question 4

DGH4 = GH * Q3.D;
Q4.Ku = margin(DGH4);

K4 = Q4.Ku / 2;

Q4.X = (K4 * DGH4)/(1 + K4*DGH4);

%step(Q4.X)

%% Question 5

Q5.p = -2*CF;
Q5.z1 = -(s+2*D+2*D*1i);
Q5.z2 = -(s+2*D-2*D*1i);

Q5.Kp = (1/Q5.p)-(Q5.z1+Q5.z2)/(Q5.z1*Q5.z2);
Q5.Kd = 1/(Q5.z1*Q5.z2) + (Q5.Kp/Q5.p);  % derivittave gain
Q5.D = (Q5.Kp + Ki*(1/s)) + Q5.Kd*(-Q5.p*s)/(s-Q5.p);

DGH5 = Q5.D * GH;

Q5.Ku = margin(DGH5);
K5 = Q5.Ku / 2;

Q4.X = (K5 * EMS)/(1 + K5*EMS);

%% Question 6

%want OS <= 10% and Ts as small as possible 

K6 = K5;

factor = 0.59; %%%%change this

Q6.p = Q5.p;
Q6.z1 = Q5.z1*factor;
Q6.z2 = Q5.z2*factor;

Q6.Kp = (1/Q6.p) - (Q6.z1+Q6.z2)/(Q6.z1*Q6.z2);
Q6.Kd = 1/(Q6.z1*Q6.z2)+(Q6.Kp/Q6.p);
Q6.D = Q6.Kp + Ki*(1/s) + Q6.Kd*(-Q6.p*s)/(s-Q6.p);

DGH6 = Q6.D *GH;
Q6.X = (K6 * DGH6)/ (1 + K6*DGH6);

stepinfo(Q6.X)