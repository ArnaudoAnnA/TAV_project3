clear all
clc
%% params
m = 1575;
Jpsi = 2875;
Jz = Jpsi;
lf = 1.3;
a = lf;
lr = 1.5;
b = lr;
Cf = 2*60e3;
Cr = 2*57e3;
l = lf+lr;

delta_max = 25;

Vxb = 50/3.6; % [m/s]


%% MATRICES

A = [0 1 0 0 
    0 -(Cf+Cr)/(m*Vxb) (Cf+Cr)/m (Cr*lr-Cf*lf)/(m*Vxb)
    0 0 0 1
    0 (Cr*lr-Cf*lf)/(Jpsi*Vxb) (-Cr*lr+Cf*lf)/Jpsi -(Cr*lr^2+Cf*lf^2)/(Jpsi*Vxb)];

B1 = [0
    Cf/m
    0
    Cf*lf/Jpsi];

B2 = [0
    (Cr*lr-Cf*lf)/(m*Vxb)-Vxb
    0
    -(Cr*lr^2+Cf*lf^2)/(Jpsi*Vxb)];

C = eye(4);

%% FB control

K = place(A, B1, [-0.10,-0.15,-0.20,-0.25]');
% K = place(A, B1, [-100,-150,-200,-250]');

Kff = m*Vxb^2/l*(lr/Cf-lf/Cr+lf/Cr*K(3))+l-lr*K(3);


%% Curvature profile
% matlab function in the simulink model
% 
% kl = zeros(2000,1);
% index = 0;
% 
% for t=0:0.1:200
%     index = int16(t*10+1);
%     if t<=25 || t>=180 || (t>=75 && t<=130)
%         kl(index) = 0;
%     elseif (t>=35 && t<=65)
%         kl(index) = -2000;
%     elseif (t>=140 && t<=170)
%         kl(index) = 2000;
%     elseif t>25 && t<35
%         kl(index) = -(t-25)*200;
%     elseif t>65 && t<75
%         kl(index) = -2000+(t-65)*200;
%     elseif t>130 && t<140
%         kl(index) = (t-130)*200;
%     elseif t>170 && t<180
%         kl(index) = 2000-(t-170)*200;
%     end
% end 
% 
% plot(kl);
%%
%open("model.slx")
%sim("model.slx")

%%
% TRANSFER FUNCTION
Ybeta = -(Cf+Cr);
Yr = 1/Vxb*(-a*Cf+b*Cr);
Ydelta = Cf;

Nbeta = (-a*Cf+b*Cr);
Nr = -1/Vxb*(a^2*Cf+b^2*Cr);
Ndelta = a*Cf;

s = tf('s');

RTF = (Ndelta*s+(Nbeta*Ydelta-Ybeta*Ndelta)/(m*Vxb))/...
    (Jz*s^2+(-Nr-Ybeta*Jz/m/Vxb)*s+(Nbeta+(Ybeta*Nr-Yr*Nbeta)/m/Vxb));

BetaTF = (Jz*Ydelta/m/Vxb*s+(-Nr*Ydelta-Ndelta*(m*Vxb-Yr))/m/Vxb)/...
    (Jz*s^2+(-Nr-Ybeta*Jz/m/Vxb)*s+(Nbeta+(Ybeta*Nr-Yr*Nbeta)/m/Vxb));

figure, step(RTF)
figure, step(BetaTF)
