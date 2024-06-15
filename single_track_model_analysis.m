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
