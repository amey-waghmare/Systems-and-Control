%% Assignment 4
% Amey Samrat Waghmare
% 203230013

% Part 1, Simulation of Nonlinear Model
clear;
clc;

SetGraphics;

load System1_Parameters;

global sys;

% Simulation Time related parameters
Simulation_Time = 100;
Samp_T = 0.1;
Ns = Simulation_Time/Samp_T;
randn('seed',0);

% Dimensions of the Dynamic System
n_st = 3;
n_op = 2;
n_ip = 2;
n_d = 1;

% Initialize Empty Matrices (*,Ns) Columns being kth Sampling instants
Xk = zeros(n_st,Ns);
Uk = zeros(n_ip,Ns);
Dk = zeros(n_d,Ns);
Yk = zeros(n_op,Ns);

% Inputs 
N_RBS = [25 30]';                 % Period of RBS
uk_ampl = [10 4]';                % Amplitude of RBS around corresponding Us
R_mat = diag([0.2 0.25].^2);  % Measurement noise Covariance

% Initial Conditions
Xk(:,1) = sys.Xs;
Yk(:,1) = sys.C_mat*Xk(:,1);

kT = zeros(1,Ns);                  % Time Array
deluk = zeros(n_ip,1);             % Corresponding to Uk
vk = zeros(n_op,1);                % Corresponding to Dk Gaussian(0,R_mat)

for k = 1:Ns-1
    kT(k) = (k-1)*Samp_T;          % Actual Time
  
    % For Random Input Signal, around Us and + or - by uk_ampl (randn ly)
    for i = 1:n_ip
        if (rem(k,N_RBS(i)) == 0)
           deluk(i) = uk_ampl(i)*sign(randn);
        end
    end
    
    Uk(:,k) = sys.Us + deluk;
    dk = normrnd(0,0.015);
    Dk(k) = sys.Ds + dk;
    
    
    % Actual Simulation Starts here
    sys.Uk = Uk(:,k);
    sys.Dk = Dk(k);
    [T,Xt] = ode45('System1_Dynamics', [0 Samp_T], Xk(:,k));       % RK solver
    Xk(:,k+1) = Xt(end,:)';
    vk = (mvnrnd(zeros(n_op,1),R_mat,1))';                         % Measurement Noise
    Yk(:,k+1) = sys.C_mat*Xk(:,k+1) + vk;
end

% At final Time Instant
kT(Ns) = Ns*Samp_T;
Uk(:,Ns) = sys.Us + deluk;
dk = sys.dk_sigma * randn;
Dk(Ns) = sys.Ds + dk;


% Figures subplot(m,n,p) Grid(m,n) and pth loc

%figure(1)
%subplot(3,1,1),plot(kT,Xk(1,:),'b-'),grid,ylabel("X1"),title('State Variables')
%subplot(3,1,2),plot(kT,Xk(2,:),'b-'),grid,ylabel("X2")
%subplot(3,1,3),plot(kT,Xk(3,:),'b-'),grid,ylabel("X3"),xlabel("Time (Min)")

%figure(2)
%subplot(2,1,1),plot(kT,Yk(1,:),'-'),grid,ylabel('Y1'),title('Measured Outputs')
%subplot(2,1,2),plot(kT,Yk(2,:),'-'),grid,ylabel('Y2'),xlabel('Time (Min)')

%figure(3)
%subplot(3,1,1),stairs(kT,Uk(1,:),'-','LineWidth',2),grid,ylabel('U1'),title('Manupulated Inputs')
%subplot(3,1,2),stairs(kT,Uk(2,:),'-','LineWidth',2),grid,ylabel('U2')
%subplot(3,1,3),stairs(kT,Dk,'-','LineWidth',2),grid,ylabel('D'),xlabel('Time (Min)')

%% Part 2: Estimation of ARX Model

y_sml = (Yk(:,1:750) - sys.Ys);
u_sml = (Uk(:,1:750) - sys.Us);

omega1_est = zeros(750-2,9);
omega2_est = zeros(750-2,9);
Y1_est = zeros(750-2,1);
Y2_est = zeros(750-2,1);

for i = 4:750
    
    omega1_est(i-3,:) = [-y_sml(1,i-3) -y_sml(1,i-2) -y_sml(1,i-1) u_sml(1,i-3) u_sml(1,i-2) u_sml(1,i-1) u_sml(2,i-3) u_sml(2,i-2) u_sml(2,i-1)];
    
    omega2_est(i-3,:) = [-y_sml(2,i-3) -y_sml(2,i-2) -y_sml(2,i-1) u_sml(1,i-3) u_sml(1,i-2) u_sml(1,i-1) u_sml(2,i-3) u_sml(2,i-2) u_sml(2,i-1)];
    
    Y1_est(i-3,:) = y_sml(1,i);
    Y2_est(i-3,:) = y_sml(2,i);
    
end

theta1_est = inv(omega1_est'*omega1_est)*omega1_est'*Y1_est;
theta2_est = inv(omega2_est'*omega2_est)*omega2_est'*Y2_est;

alpha_1est = theta1_est(1:3);
beta1_1est = theta1_est(4:6);
beta2_1est = theta1_est(7:9);

alpha_2est = theta2_est(1:3);
beta1_2est = theta2_est(4:6);
beta2_2est = theta2_est(7:9);

% Now the Models are
phi_1_est = zeros(3,3);
phi_1_est(:,1) = -1*flip(alpha_1est);
phi_1_est(:,2) = [1;0;0]; phi_1_est(:,3) = [0;1;0];

L_1_est = -1*flip(alpha_1est);

gamma_1_est = zeros(3,2);
gamma_1_est(:,1) = flip(beta1_1est);
gamma_1_est(:,2) = flip(beta2_1est);


phi_2_est = zeros(3,3);
phi_2_est(:,1) = -1*flip(alpha_2est);
phi_2_est(:,2) = [1;0;0]; phi_2_est(:,3) = [0;1;0];

L_2_est = -1*flip(alpha_2est);

gamma_2_est = zeros(3,2);
gamma_2_est(:,1) = flip(beta1_2est);
gamma_2_est(:,2) = flip(beta2_2est);

% Combined Estimated model
phi_est = [phi_1_est zeros(3,3);zeros(3,3) phi_2_est];
gamma_est = [gamma_1_est;gamma_2_est];
L_est = [L_1_est zeros(size(L_2_est));zeros(size(L_1_est)) L_2_est];
C_est = [1 0 0 0 0 0;0 0 0 1 0 0];

%% Part 3, Validation (Prediction,Simulation)
y_val = (Yk(:,751:end) - sys.Ys);
u_val = (Uk(:,751:end) - sys.Us);
n_state_est = length(phi_est);
xk_pred = zeros(n_state_est, 251);
xk_sim = zeros(n_state_est, 251);
yk_pred = zeros(n_op, 251);
yk_sim = zeros(n_op, 251);

for k=1:250
    ek = y_val(:,k) - (C_est*xk_pred(:,k));
    xk_pred(:,k+1) = phi_est*xk_pred(:,k) + gamma_est*u_val(:,k) + L_est*ek;
    yk_pred(:,k+1) = (C_est*xk_pred(:,k+1));
    
    xk_sim(:,k+1) = phi_est*xk_sim(:,k) + gamma_est*u_val(:,k);
    yk_sim(:,k+1) = C_est * xk_sim(:,k+1);
end
%% Part 4, Input and Output plots

Kt = 75:Samp_T:100;
figure(1)
subplot(2,1,1),plot(Kt,Yk(1,750:end)-sys.Ys(1),':r',Kt,yk_pred(1,:),'b',Kt,yk_sim(1,:),'g'),grid,ylabel('yk, yk pred, yk sim'),title('Measured Outputs')
subplot(2,1,2),plot(Kt,Yk(2,750:end)-sys.Ys(2),':r',Kt,yk_pred(2,:),'b',Kt,yk_sim(2,:),'g'),grid,ylabel('yk, yk pred, yk sim'),xlabel('Time (Min)')

figure(2)
subplot(2,1,1),plot(Kt(1:end-1),u_val(1,:)),grid,ylabel('u_1'),title('Inputs')
subplot(2,1,2),plot(Kt(1:end-1),u_val(2,:)),grid,ylabel('u_2'),xlabel('Time (min)')
%% Part 5, Comparision of Step Response

load System1_ContLinmod;

continous_model_original = ss(A_mat, B_mat, C_mat, D_mat);
discrete_model_original = c2d(continous_model_original, Samp_T);

discrete_model_estimated = ss(phi_est, gamma_est, C_est, zeros(2,2),Samp_T);

figure(3)
step(discrete_model_original,'r--', discrete_model_estimated),legend({'Original Model','Estimated Model'})
