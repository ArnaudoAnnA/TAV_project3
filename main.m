clear all
clc
%% Parameters
m = 1575;
Jpsi = 2875;
Jz = Jpsi;
lf = 1.3;
a = lf;
lr = 1.5;
b = lr;
Cf = 2*60e3;
Cr = 2*57e3;% under
%Cr = 57e3; % over
%Cr = Cf*a/b; % neutral
l = lf+lr;
CR=Cr;
CF=Cf;
L=l;
tau_s=1;

delta_max = 25;

Vxb = 80/3.6; % [m/s]

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

%% SINGLE TRACK model analysis
% Poles, Damping Factors, Natural Frequencies versus velocity
% vehicle speed
vel=[[1:0.1:20],[20:0.5:50],[50:5:450]]/3.6;
delta_f = 1*pi/180;
delta_f = 1;
delta_r = 0;
% matrix initialization 
POLES=[zeros(length(vel),2)];       % Poles
ZETA=[zeros(length(vel),2)];        % damping fators
FREQ=[zeros(length(vel),2)];        % natural frequencies
DET_A=[zeros(length(vel),1)];       % determinant of A
TR_A=[zeros(length(vel),1)];        % trace of A
X_r = [zeros(length(vel),2)];
Y_r = [zeros(length(vel),6)];
ay_r = [zeros(length(vel),1)];
u_r = [delta_f;0];
values = [zeros(length(vel),2)];

% Definition of state space variables
StateNames ={'\beta','r'};
InputNames={'\delta_F','\delta_R'};
OutputNames={'\beta','r','\rho','\alpha_F','\alpha_R','a_y'};

% for cycle (vehicle speed)
for k=1:length(vel)
    Vv=vel(k);      % vehicle speed
    % state space matrices: A,B,C,D
    A=[(-CF-CR)/(m*Vv),(-CF*a+CR*b-m*Vv^2)/(m*Vv^2);
        (-CF*a+CR*b)/Jz,(-CF*a^2-CR*b^2)/(Jz*Vv)];
    B=[CF/(m*Vv) CR/(m*Vv);
        (CF*a/Jz) -(CR*b/Jz)];
    C= [1,0
        0,1
        (-CR-CF)/(m*Vv^2),(-CF*a+CR*b)/(m*Vv^3)
        -1, -a/Vv
        -1, b/Vv
        (-CR-CF)/(m),(-CF*a+CR*b)/(m*Vv)];
    D = [0 0;
        0 0;
        CF/(m*Vv^2) CR/(m*Vv^2)
        1 0
        0 1
        CF/m CR/m];
    % state space system
    G=ss(A,B,C,D);
    % Wn, Z e P: fnatural frequencies, damping factors and poles
    [Wn,Z,P]=damp(G);
    % storing values in matrices (k-th column)
    POLES(k,:)=P;
    ZETA(k,:)=Z;
    FREQ(k,:)=Wn;
    % Determinant and trace of matrix A
    DET_A(k)=det(A);
    TR_A(k)=trace(A);
    % steady state response to delta_F
    X_r(k,:) = -A^-1*B*u_r;
    Y_r(k,:) = C*X_r(k,:)'+D*u_r;
    ay_r(k,:) = X_r(k,2)'*Vv;
    [V,D] = eig(A);
    values(k,:)=eig(A);
end

% poles
P1=POLES(:,1);
P2=POLES(:,2);
% real and imaginary part (first pole)
P1_Real=real(P1);
P1_Im=imag(P1);
% real and imaginary part (second pole)
P2_Real=real(P2);
P2_Im=imag(P2);

fig_1=figure('Name','Poles, Damping, Natural Frequencies','NumberTitle','off','PaperType','A4');
figure(fig_1)
subplot(2,2,[3 4])
% Plot eigenvalues vs speed
scatter(P1_Real,P1_Im,[],[1:1:length(vel)])
hold on
grid on
scatter(P2_Real,P2_Im,[],[1:1:length(vel)])
xlabel('Real'),ylabel('Im'),title('Poles'),
% xlim([-200 20]),ylim([-10 10])
plot([0,0],ylim,'--k')
hold on
colorbar %axis equal
set(gca,'FontName','Times New Roman','FontSize',14)

subplot(2,2,1),plot(vel*3.6,FREQ,'LineWidth',2)
title('Natural Frequencies'),xlabel('vel [km/h]'),ylabel('Wn [rad/s]')
grid on
set(gca,'FontName','Times New Roman','FontSize',14)

subplot(2,2,2),plot(vel*3.6,ZETA,'square')
xlabel('vel [km/h]'),ylabel('Damping Ratio'),title('Damping Ratio');
grid on
%plot(POLES,'x'),hold on,
set(gca,'FontName','Times New Roman','FontSize',14)

% Plot determinant versus speed
fig_11=figure('Name','Determinant of [A]','NumberTitle','off','PaperType','A4');
figure(fig_11)
hold all
plot(vel*3.6,DET_A,'linewidth',2),
plot(vel*3.6,zeros(length(vel),1),'--k')
xlabel('vel [km/h]'),ylabel('det(A)'),title('$\lambda_1 \lambda_2$','interpreter','latex','Fontsize',18)
ylim([-10 100])
set(gca,'FontName','Times New Roman','FontSize',14)
figure, 
plot(vel*3.6,values)
figure('Name','steady state response vs. velocity')
subplot(2,2,1)
plot(vel*3.6,Y_r(:,1),'o','linewidth',2,'Displayname','$\beta$'); %xlabel('vel [km/h]'),%ylabel('Y_r'),
title('$\beta /\delta_F$','interpreter','latex','Fontsize',18); hold on; 
plot(xlim,[0 0],'--k'); grid on
% legend('\beta [rad]','r [rad/s]','\rho [1/m]','\alpha_f [rad]','\alpha_r [rad]')

subplot(2,2,2)
plot(vel*3.6,Y_r(:,2),'o','linewidth',2)
title('$r /\delta_F$','interpreter','latex','Fontsize',18); grid on

subplot(2,2,3)
plot(vel*3.6,Y_r(:,3),'o','linewidth',2); title('$\rho /\delta_F$','interpreter','latex','Fontsize',18); xlabel('vel [km/h]'); grid on

subplot(2,2,4)
plot(vel*3.6,Y_r(:,4:5),'o','linewidth',1);legend('$\alpha_f /\delta_F$','$\alpha_r /\delta_F$','Location','best','interpreter','latex','Fontsize',18)
xlabel('vel [km/h]')
grid on

figure('name','lateral acceleration')
plot(vel*3.6,Y_r(:,6),'o','linewidth',2); 
title('$a_y /\delta_F$','interpreter','latex','Fontsize',18)
xlabel('vel [km/h]'); grid on; xlabel('vel [km/h]')
delta_meno_delta_0 = u_r(1)-L*Y_r(:,3);
beta_meno_beta_0 = Y_r(:,1)-b*Y_r(:,3);

ay_r = Y_r(:,6);
figure('name','understeering diagram')
subplot(1,2,1)
plot(ay_r,delta_meno_delta_0*180/pi*tau_s,'linewidth',2)
xlabel('a_y [m/s^2]'); grid on; 
ylabel('$\delta-\delta_0$ [deg]','interpreter','latex','Fontsize',18); xlim([0 5])
subplot(1,2,2)
plot(ay_r,beta_meno_beta_0*180/pi,'linewidth',2); 
xlabel('a_y [m/s^2]'); grid on; 
ylabel('$\beta-\beta_0$ [deg]','interpreter','latex','Fontsize',18); xlim([0 5])

K_us_num = (delta_meno_delta_0(end)- delta_meno_delta_0(1))/(ay_r(end)-ay_r(1))
K_beta_num = (beta_meno_beta_0(end)- beta_meno_beta_0(1))/(ay_r(end)-ay_r(1))

pause

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
K = acker(A, B1, [-2, -2, -2, -2]');
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
    %gain_matrix(k,:) = acker(A, B1, [-5, -5, -5, -5]');
    sys=ss(Ak,B1k,Ck,zeros(4,1));
    gain_matrix(k,:) = lqr(sys,Q,R);

    Kff_matrix(k) = m*Vxb^2/l*(lr/Cf-lf/Cr+lf/Cr*gain_matrix(k,3))+l-lr*gain_matrix(k,3);
end

figure, plot(vel, gain_matrix)
legend('k(1)', 'k(2)', 'k(3)', 'k(4)')
figure, plot(vel, gain_matrix(:,1))
figure, plot(vel, gain_matrix(:,2))
figure, plot(vel, gain_matrix(:,3))
figure, plot(vel, gain_matrix(:,4))
figure, plot(vel, Kff_matrix);


%% State space matrices for simulation
% slide 29 TAV_AS_2.6
CF=Cf;
CR=Cr;
Vv=Vxb;
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
t_end_sim=200;


%% SIMULATION
open("model.slx")
out=sim("model.slx");
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
OutputNames={'\beta','r','\rho','\alpha_F','\alpha_R','a_y'};
L=l;
tau_s=1;
figure('Name','steering angle')
hold all; grid on
plot(out.delta,'LineWidth',2),xlabel('time [s]')
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
dt_sim = 1e-3
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


