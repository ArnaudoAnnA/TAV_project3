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

Vxb = 80/3.6; % [m/s]

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

% K = place(A, B1, [-0.10,-0.15,-0.20,-0.25]');
% K = place(A, B1, [-100,-150,-200,-250]');
K = acker(A, B1, [-1, -1, -1, -1]');

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
Vx_50=50/3.6;
Vx_120=120/3.6;
Ybeta = -(Cf+Cr);
Yr = 1/Vxb*(-a*Cf+b*Cr);
Yr_50 = 1/Vx_50*(-a*Cf+b*Cr);
Yr_120 = 1/Vx_120*(-a*Cf+b*Cr);
Ydelta = Cf;


Nbeta = (-a*Cf+b*Cr);
Nr = -1/Vxb*(a^2*Cf+b^2*Cr);
Ndelta = a*Cf;

Nbeta_50 = (-a*Cf+b*Cr);
Nr_50 = -1/Vx_50*(a^2*Cf+b^2*Cr);
Ndelta_50 = a*Cf;

Nbeta_120 = (-a*Cf+b*Cr);
Nr_120 = -1/Vx_120*(a^2*Cf+b^2*Cr);
Ndelta_120 = a*Cf;

s = tf('s');

RTF_50 = (Ndelta_50*s+(Nbeta_50*Ydelta-Ybeta*Ndelta_50)/(m*Vx_50))/...
    (Jz*s^2+(-Nr-Ybeta*Jz/m/Vx_50)*s+(Nbeta_50+(Ybeta*Nr_50-Yr*Nbeta_50)/m/Vx_50));

BetaTF_50 = (Jz*Ydelta/m/Vx_50*s+(-Nr*Ydelta-Ndelta_50*(m*Vx_50-Yr))/m/Vx_50)/...
    (Jz*s^2+(-Nr_120-Ybeta*Jz/m/Vx_50)*s+(Nbeta_50+(Ybeta*Nr_50-Yr_50*Nbeta_50)/m/Vx_50));

RTF_120 = (Ndelta_120*s+(Nbeta_120*Ydelta-Ybeta*Ndelta_120)/(m*Vx_120))/...
    (Jz*s^2+(-Nr-Ybeta*Jz/m/Vx_120)*s+(Nbeta_120+(Ybeta*Nr_120-Yr*Nbeta_120)/m/Vx_120));

BetaTF_120 = (Jz*Ydelta/m/Vx_120*s+(-Nr_120*Ydelta-Ndelta_120*(m*Vx_120-Yr))/m/Vx_120)/...
    (Jz*s^2+(-Nr_120-Ybeta*Jz/m/Vx_120)*s+(Nbeta_120+(Ybeta*Nr_120-Yr_120*Nbeta_120)/m/Vx_120));

figure, step(RTF_50)
hold on
step(RTF_120)

figure, step(BetaTF_50)
hold on
step(BetaTF_120)


%% State space matrices for simulation
CF=Cf;
CR=Cr;
Vv=Vxb;
A_sim=[(-CF-CR)/(m*Vv),(-CF*a+CR*b-m*Vv^2)/(m*Vv^2);
    (-CF*a+CR*b)/Jz,(-CF*a^2-CR*b^2)/(Jz*Vv)];
B_sim=[CF/(m*Vv) CR/(m*Vv);
    (CF*a/Jz) -(CR*b/Jz)];
C_sim = [1,0
    0,1
    (-CR-CF)/(m*Vv^2),(-CF*a+CR*b)/(m*Vv^3)
    -1, -a/Vv
    -1, b/Vv
    (-CR-CF)/(m),(-CF*a+CR*b)/(m*Vv)];
D_sim = [0 0;
    0 0;
    CF/(m*Vv^2) CR/(m*Vv^2)
    1 0
    0 1
    CF/m CR/m];
t_end_sim=200;


pause


%% Bode plot 
% clear A B C D
%     A=[(-CF-CR)/(m*Vv),(-CF*a+CR*b-m*Vv^2)/(m*Vv^2);
%         (-CF*a+CR*b)/Jz,(-CF*a^2-CR*b^2)/(Jz*Vv)];
%     B=[CF/(m*Vv) CR/(m*Vv);
%         (CF*a/Jz) -(CR*b/Jz)];
%     C = [1,0
%         0,1
%         (-CR-CF)/(m*Vv^2),(-CF*a+CR*b)/(m*Vv^3)
%         -1, -a/Vv
%         -1, b/Vv
%         (-CR-CF)/(m),(-CF*a+CR*b)/(m*Vv)];
%     D = [0 0;
%         0 0;
%         CF/(m*Vv^2) CR/(m*Vv^2)
%         1 0
%         0 1
%         CF/m CR/m];
% sys=ss(A,B,C,D);
% fig_3=figure('Name','Bode Plots','NumberTitle','off','PaperType','A4');
% figure(fig_3)
% 
% figure; step(sys); figure; ltiview(sys)
% h = bodeplot(sys);
% % % Change units to Hz and make phase plot invisible
% setoptions(h,'FreqUnits','Hz','PhaseVisible','on','Grid','On','Xlim',[0.1 10],'MagScale','linear','MagUnits','abs'); %,'PhaseWrapping','on'
% set(findall(gcf,'-property','FontSize'),'FontSize',13)


%% POST PROCESSING
L=l;
tau_s=1;
figure('Name','steering angle')
hold all; grid on
plot(out.delta,'LineWidth',2),xlabel('time [s]'),
%hold on; plot(L*out.ro*180/pi*tau_s,'--k'); 
%legend('\delta','\delta_0','Fontsize',18,'location','best')
%title('Steering Angle \delta_s')

%--------- Plot beta and psi_dot
figure('Name','States')
hold all; grid on
subplot(2,1,1),plot(out.beta,'LineWidth',2),xlabel('time [s]')
hold on
plot(b*out.ro*180/pi,'--k'); 
legend('\beta','\beta_0','Fontsize',16,'location','best')
title('slip angle \beta [deg]','Fontsize',16)
grid on
subplot(2,1,2)
plot(out.r,'LineWidth',2)
hold on
plot(out.delta/180*pi,'LineWidth',2),xlabel('time [s]'),
title('r [deg/s]','Fontsize',16), xlabel('time [s]')
ylabel('')
legend('r [deg/s]','\delta [rad]','Fontsize',16,'location','best')
%ylim([-150 150])
grid on

% --------- Plot ay vs t
figure('Name','a_y(t)')
hold all; grid on
plot(out.ay,'LineWidth',2)
title('Lateral Acceleration a_y [m/s^2]','Fontsize',16)
xlabel('time [s]')
ylabel('')
% ---- plot beta beta_dot
figure('Name','\beta, \beta_dot')
hold all; grid on
plot(out.beta.Data,out.beta_dot(:,2),'LineWidth',2) %title('Lateral Acceleration a_y [m/s^2]'),
xlabel('\beta [deg]','Fontsize',18)
ylabel('\beta_{dot} [rad/s]','Fontsize',18)


% %% --------- Plot beta vs ay
% figure('Name','beta vs ay')
% % plot(a_y(:,2),Beta(:,2))
% scatter(out.ay(:,2),out.Beta(:,2),[],out.ay(:,1)); colorbar
% xlabel('a_y [m/s^2]')
% ylabel('\beta [deg]'); grid on
% text(11,2.2,['time[s]'])
% set(gca,'FontName','Times New Roman','FontSize',16)
% %--------- Plot delta vs ay
% figure('Name','delta vs ay')
% % plot(a_y(:,2),Beta(:,2))
% scatter(out.ay(:,2),out.delta_rad(:,2)*180/pi,[],out.ay(:,1)); colorbar
% xlabel('a_y [m/s^2]'); ylabel('\delta_{vol} [deg]'); grid on
% text(11,22,['time[s]'])
% set(gca,'FontName','Times New Roman','FontSize',16)
% %--------- Plot delta-delta0 vs ay
% figure('Name','delta-delta_0 vs ay')
% delta0 =L*ro.Data;
% % plot(a_y(:,2),Beta(:,2))
% plot(out.ay(:,2),(out.delta_rad(:,2)-delta0*tau_s)*180/pi,'linewidth',2); 
% xlabel('a_y [m/s^2]'); ylabel('\delta_{vol}-\delta_0 [deg]'); grid on
% set(gca,'FontName','Times New Roman','FontSize',16)


% --------- Plot alpha_F e alpha_R 
figure('Name','alphaF e R'); hold all
plot(out.alfaF*180/pi,'LineWidth',2); plot(out.alfaR*180/pi,'LineWidth',2); 
xlabel('time [s]'); ylabel('\alpha [deg]'); grid on; 
legend('\alpha_F','\alpha_R','Fontsize',16,'location','best')
set(gca,'FontName','Times New Roman','FontSize',14)
legend({},'FontSize',16)

%--------- Plot curvature
figure('Name','rho'); hold all
plot(out.ro,'LineWidth',2); 
xlabel('time [s]'); ylabel('\rho [1/m]'); grid on; 
set(gca,'FontName','Times New Roman','FontSize',14)
%--------- Plot trajectory
spost_x=out.Var_trajectory(:,1);
spost_y=out.Var_trajectory(:,2);

% figure
% hold all; grid on
% % plot(spost_x,spost_y,'LineWidth',2)
% scatter(spost_x,spost_y,[],a_y(:,1))
% title('trajectory'),axis equal,xlabel('X [m]'),ylabel('Y[m]');colorbar
% text(49,12,['time[s]'])

