clear all
clc

load_params;

% check under/over steering
if CR*b-CF*a>0
    disp('============ understeering ============')
elseif CF*a-CR*b==0
    disp('============ neutral vehicle ============')
else
    disp('============ oversteering ============')
    V_cr = sqrt(CF*CR*L^2/(m*(a*CF-b*CR)))*3.6;
    disp(['critical speed: ',num2str(round(V_cr*10)/10),' km/h'])
end

%single_track_model_analysis;

%% TRANSFER FUNCTION
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
% hold on
step(RTF_120)

figure, step(BetaTF_50)
% hold on
step(BetaTF_120)

close all
%% MATRICES
% used in the control model
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


%% FB and FF controllers
% Matrix K with pole placement
% K = place(A, B1, [-0.10,-0.15,-0.20,-0.25]');
% K = place(A, B1, [-100,-150,-200,-250]');
K = acker(A, B1, [-10, -10, -10, -10]');
Kff = m*Vxb^2/l*(lr/Cf-lf/Cr+lf/Cr*K(3))+l-lr*K(3);


% Matrix K with linear quadratic set-up
Q = [1,0,0,0;...
    0,1,0,0;...
    0,0,1,0;...
    0,0,0,1];
R = 1;
D=zeros(4,1);
sys=ss(A,B1,C,D);
K1=lqr(sys,Q,R);
Kff1 = m*Vxb^2/l*(lr/Cf-lf/Cr+lf/Cr*K1(3))+l-lr*K1(3);

%% Plot how K matrix varies wrt vehicle speed
vel=linspace(1,90)/3.6;
gain_matrix = zeros(length(vel), 4);
Kff_matrix = zeros(length(vel));
for k=1:length(vel)
    Vx=vel(k);      % vehicle speed
    % state space matrices: A,B,C,D
    Ak = [0 1 0 0 
    0 -(Cf+Cr)/(m*Vx) (Cf+Cr)/m (Cr*lr-Cf*lf)/(m*Vx)
    0 0 0 1
    0 (Cr*lr-Cf*lf)/(Jpsi*Vx) (-Cr*lr+Cf*lf)/Jpsi -(Cr*lr^2+Cf*lf^2)/(Jpsi*Vx)];

    B1k = [0
        Cf/m
        0
        Cf*lf/Jpsi];
    
    B2k = [0
        (Cr*lr-Cf*lf)/(m*Vx)-Vx
        0
        -(Cr*lr^2+Cf*lf^2)/(Jpsi*Vx)];
    
    Ck = eye(4);
    %gain_matrix(k,:) = acker(Ak, B1k, [-10, -10, -10, -10]');
    sys=ss(Ak,B1k,Ck,zeros(4,1));
    gain_matrix(k,:) = lqr(sys,Q,R);

    Kff_matrix(k) = m*Vxb^2/l*(lr/Cf-lf/Cr+lf/Cr*gain_matrix(k,3))+l-lr*gain_matrix(k,3);
end

figure, plot(vel, gain_matrix)
legend('k(1)', 'k(2)', 'k(3)', 'k(4)')
figure, plot(vel, gain_matrix(:,1)), xlabel('V(m/s)'), title('K(1)')
figure, plot(vel, gain_matrix(:,2)), xlabel('V(m/s)'), title('K(2)')
figure, plot(vel, gain_matrix(:,3)), xlabel('V(m/s)'), title('K(3)')
figure, plot(vel, gain_matrix(:,4)), xlabel('V(m/s)'), title('K(4)')
figure, plot(vel, Kff_matrix);

%% State space matrices for simulation
A_sim=[(-CF-CR)/(m*Vv),(-CF*a+CR*b-m*Vv^2)/(m*Vv^2);
    (-CF*a+CR*b)/Jz,(-CF*a^2-CR*b^2)/(Jz*Vv)];
B_sim=[CF/(m*Vv) CR/(m*Vv); %aggiunto -
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


%% SIMULATION
open("model.slx")
out=sim("model.slx");

%% POST PROCESSING
OutputNames={'\beta','r','\rho','\alpha_F','\alpha_R','a_y'};

figure('Name','steering angle')
hold all; grid on
plot(out.delta,'LineWidth',2),xlabel('time [s]')
%hold on; plot(L*out.ro*180/pi*tau_s,'--k'); 
legend('\delta','\delta_0','Fontsize',18,'location','best')
title('Steering Angle \delta_s')

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

%% --------- Plot alpha_F e alpha_R 
figure('Name','alphaF e R'); hold all
plot(out.alfaF*180/pi,'LineWidth',2); plot(out.alfaR*180/pi,'LineWidth',2); 
xlabel('time [s]'); ylabel('\alpha [deg]'); grid on; 
legend('\alpha_F','\alpha_R','Fontsize',16,'location','best')
set(gca,'FontName','Times New Roman','FontSize',14)
legend({},'FontSize',16)

%% --------- Plot curvature
figure('Name','rho'); hold all
plot(out.ro,'LineWidth',2); 
xlabel('time [s]'); ylabel('\rho [1/m]'); grid on; 
set(gca,'FontName','Times New Roman','FontSize',14)
%--------- Plot trajectory
spost_x=out.Var_trajectory(:,1);
spost_y=out.Var_trajectory(:,2);

figure
hold all; grid on
% plot(spost_x,spost_y,'LineWidth',2)
scatter(spost_x,spost_y,[],out.ay.Data(:,1))
%title('trajectory'),xlabel('X [m]'),ylabel('Y[m]');colorbar
title('trajectory'),axis equal,xlabel('X [m]'),ylabel('Y[m]');colorbar

text(49,12,['time[s]'])

%% Vehicle Trajectory and Animation
F_Size=13;

figure('Name','Vehicle CG location','NumberTitle','off','PaperType','A4');
hold all; grid on
plot(out.Var_trajectory(:,1),out.Var_trajectory(:,2),'Linewidth',2);
title('CG Trajectory'); axis equal
set(gca,'FontName','Times New Roman','FontSize',F_Size)
xlabel('X [m]'); ylabel('Y [m]')
X_G = out.Var_trajectory(:,1);
Y_G = out.Var_trajectory(:,2);

% figure
%hold on
axis equal
dt_sim = 1e-3;
dt_frame = 0.1; % [s] vehicle frame refresh (time) interval 
decim_frame = dt_frame/dt_sim;

% cycle to show vehicle motion
for cont1=1:decim_frame:length(out.Psi)
    X = X_G(cont1)+[-b a a -b]';
    Y = Y_G(cont1)+[-1 -1 1 1]';
    vert = [X,Y];
    fac = [1 2 3 4];
    hVeicolo = patch('Faces',fac,'Vertices',vert,'FaceColor','red','FaceAlpha',.5);
    direction = [0 0 1];
    xlim([X(1)-5 X(2)+3])
    ylim([Y(1)-5 Y(3)+3])
    
    x0 = X_G(cont1);
    y0 = Y_G(cont1);
    z0 = 0;
    ORIGIN = [x0,y0,z0];
    rotate(hVeicolo,direction,out.Psi(cont1,2),ORIGIN);
    pause(0.1)
end

plot(out.Var_trajectory(:,1),out.Var_trajectory(:,2),'Linewidth',2);
axis auto
