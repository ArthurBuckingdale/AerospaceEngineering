clc
clear all
close all

D = readtable('systemesDeNavigation.csv'); % attitude reel


% values are in meters
xm = D{:,3}*1000;
ym = D{:,4}*1000;
zm = D{:,5}*1000;
vxm = D{:,6}*1000;
vym = D{:,7}*1000;
vzm = D{:,8}*1000;

% time is in seconds
%time = datetime(D{:,9});
timeStamps = D{:,9};
for i = 1:length(D{:,9})
    time(i) = (str2num(timeStamps{i}(end-11:end-10)) + str2num(timeStamps{i}(end-8:end-7))/60 + str2num(timeStamps{i}(end-5:end))/3600) * 3600;
end
% variance
sigma_velocity = 50; %(m/s)
sigma_measure = 1; % (m/s)
sigma_position = 1000; %meters

%generate noise
W = [sigma_position 0 0 0 0 0;
    0 sigma_position 0 0 0 0;
    0 0 sigma_position 0 0 0;
    0 0 0 sigma_velocity 0 0;
    0 0 0 0 sigma_velocity 0;
    0 0 0 0 0 sigma_velocity;];
V = [1 0 0;
    0 1 0;
    0 0 1];

%matrices
G = [1 0 0 0 0 0;
    0 1 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1;];

D = [1 0 0;
    0 1 0;
    0 0 1];

P = [sigma_position^3 0 0 0 0 0;
    0 sigma_position^3 0 0 0 0;
    0 0 sigma_position^3 0 0 0;
    0 0 0 sigma_velocity^3 0 0;
    0 0 0 0 sigma_velocity^3 0;
    0 0 0 0 0 sigma_velocity^3;]; % initial guess of P

% C = [0 0 0 1 0 0;
%     0 0 0 0 1 0;
%     0 0 0 0 0 1];


stateEst = zeros(6,length(xm)); % length will probably change when we merge code
x0 = [xm(1);ym(1);zm(1);vxm(1);vym(1);vzm(1)]; % initiale value
stateEst(:,1) = x0;
stateNoCorr = stateEst; % the state without any corrections
stateNoEkf = stateEst;

trueState = [xm,ym,zm,vxm,vym,vzm]; % true measurements

interval = zeros(length(xm),1); % time interval between measurements


% EKF
Kmax = 99; % %%%%%%%%%%%_MODIFY__%%%%%%%%%%%% change this if you want to check more or less values, Max is 99;

time_vector = zeros(Kmax,1);


X_moins(1,:) = trueState(1,:);
P_moins(:,:,1) = P;
for i = 1:100
    angularVelAiding(i,:) = cross(trueState(i,1:3),trueState(i,4:6))/(norm(trueState(i,1:3)).^2);
    regVel(i,:) = cross(angularVelAiding(i,:),trueState(i,1:3));
end
angularVelAiding = angularVelAiding + 0.00001 * randn(100,3);

for k = 1:Kmax
    disp(k)
    interval(k) = time(k+1) - time(k);
    time_vector(k) = time(k) - time(1);
    
    % linearization and gain calculation for angular velocities
    C = [0 trueState(k,6) -trueState(k,5) 0 -trueState(k,3) trueState(k,2);
        -trueState(k,6) 0 trueState(k,4) trueState(k,3) 0 -trueState(k,1);
        trueState(k,5) -trueState(k,4) 0 -trueState(k,2) trueState(k,1) 0]./(norm(trueState(k,1:3)).^2);
    K = P_moins(:,:,k)*C'/(C*P_moins(:,:,k)*C' + D*V*D');
    F = findJacobian(X_moins(k,:));
    
    % measurement update step
    K*(angularVelAiding(k,:)' - C*X_moins(k,:)')
    X_plus = X_moins(k,:)' + K*(angularVelAiding(k,:)' - C*X_moins(k,:)'); % update state
    P = (eye(6) - K*C)*P_moins(:,:,k); % update covariance estimate
    
    if k == Kmax 
        break
    end
    % prediction of next vector
    [~,integratedState] = ode45(@findStateDot, [time(k) time(k+1)], X_plus);
    X_moins(k+1,:) = integratedState(end,:)'; % state prediction
    
    % prediction of covariance
    [~,integratedCov] = ode45(@(~,integratedCov) findCovarianceDot(P,F,G,W), [time(k) time(k+1)],P);
    P_moins(:,:,k+1) = reshape(integratedCov(end,:),6,6); % covariance predecition
    
    
end

% prediction of next vector
[t,integratedState] = ode45(@findStateDot, [time(1) time(end)], trueState(1,:));
figure
subplot(2,2,1)
hold on
plot(t,integratedState(:,1))
plot(time(1:Kmax),trueState(1:Kmax,1),'*')
hold off
title('Inertial Integration X and GPS')
xlabel('Time (s)')
ylabel('Position (m)')
legend('Inertial','GPS')
subplot(2,2,2)
hold on
plot(t,integratedState(:,2))
plot(time(1:Kmax),trueState(1:Kmax,2),'*')
hold off
title('Inertial Integration Y and GPS')
xlabel('Time (s)')
ylabel('Position (m)')
legend('Inertial','GPS')
subplot(2,2,3)
hold on
plot(t,integratedState(:,3))
plot(time(1:Kmax),trueState(1:Kmax,3),'*')
hold off
title('Inertial Integration Z and GPS')
xlabel('Time (s)')
ylabel('Position (m)')
legend('Inertial','GPS')


figure
subplot(2,3,1)
hold on
plot(time(1:Kmax),X_moins(:,1),'-*')
plot(time(1:Kmax),trueState(1:Kmax,1),'-*')
title('GPS and Kalman Corrected Measures in Position X')
xlabel('Time (s)')
ylabel('Position (m)')
legend({'$\hat{x}$','GPS'},'Interpreter','latex')
hold off
subplot(2,3,2)
hold on
plot(time(1:Kmax),X_moins(:,2),'-*')
plot(time(1:Kmax),trueState(1:Kmax,2),'-*')
hold off
title('GPS and Kalman Corrected Measures in Position Y')
xlabel('Time (s)')
ylabel('Position (m)')
legend({'$\hat{y}$','GPS'},'Interpreter','latex')
subplot(2,3,3)
hold on
plot(time(1:Kmax),X_moins(:,3),'-*')
plot(time(1:Kmax),trueState(1:Kmax,3),'-*')
hold off
title('GPS and Kalman Corrected Measures in Position Z')
xlabel('Time (s)')
ylabel('Position (m)')
legend({'$\hat{z}$','GPS'},'Interpreter','latex')
subplot(2,3,4)
hold on
plot(time(1:Kmax),X_moins(:,4),'-*')
plot(time(1:Kmax),trueState(1:Kmax,4),'-*')
hold off
title('GPS and Kalman Corrected Measures in Velocity Vx')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
legend({'$\hat{v}_x$','GPS'},'Interpreter','latex')
subplot(2,3,5)
hold on
plot(time(1:Kmax),X_moins(:,5),'-*')
plot(time(1:Kmax),trueState(1:Kmax,5),'-*')
hold off
title('GPS and Kalman Corrected Measures in Velocity Vy')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
legend({'$\hat{v}_y$','GPS'},'Interpreter','latex')
subplot(2,3,6)
hold on
plot(time(1:Kmax),X_moins(:,6),'-*')
plot(time(1:Kmax),trueState(1:Kmax,6),'-*')
hold off
title('GPS and Kalman Corrected Measures in Velocity Vz')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
legend({'$\hat{v}_z$','GPS'},'Interpreter','latex')

figure
subplot(2,3,1)
plot(time(1:Kmax),sqrt(reshape(P_moins(1,1,:),1,[])),'-*')
title('Covariance in X Position')
xlabel('Time (s)')
ylabel('Covariance')
subplot(2,3,2)
plot(time(1:Kmax),sqrt(reshape(P_moins(2,2,:),1,[])),'-*')
title('Covariance in Y Position')
xlabel('Time (s)')
ylabel('Covariance')
subplot(2,3,3)
plot(time(1:Kmax),sqrt(reshape(P_moins(3,3,:),1,[])),'-*')
title('Covariance in Z Position')
xlabel('Time (s)')
ylabel('Covariance')
subplot(2,3,4)
plot(time(1:Kmax),sqrt(reshape(P_moins(4,4,:),1,[])),'-*')
title('Covariance in X Velocity')
xlabel('Time (s)')
ylabel('Covariance')
subplot(2,3,5)
plot(time(1:Kmax),sqrt(reshape(P_moins(5,5,:),1,[])),'-*')
title('Covariance in Y Velocity')
xlabel('Time (s)')
ylabel('Covariance')
subplot(2,3,6)
plot(time(1:Kmax),sqrt(reshape(P_moins(6,6,:),1,[])),'-*')
title('Covariance in Z Velocity')
xlabel('Time (s)')
ylabel('Covariance')

% disp('code paused')
% pause
% 
% for k = 2:Kmax % we correct the value only when new measurement is available
%     
%     [~,integratedState] = ode45(@findStateDot, [0 interval(k)], stateEst(:,k-1));
%     X_moins = integratedState(end,:)';
%     
%     stateNoEkf(:,k) = X_moins;
%     
% end
% 
% 
% % error for the corrected state
% erreur_x = sqrt((xm(1:k) - stateEst(1,1:k)').^2); % in m
% erreur_y = sqrt((ym(1:k) - stateEst(2,1:k)').^2); % in m
% erreur_z = sqrt((zm(1:k) - stateEst(3,1:k)').^2); % in m
% erreur_vx = sqrt((vxm(1:k) - stateEst(4,1:k)').^2); % in m/s
% erreur_vy = sqrt((vym(1:k) - stateEst(5,1:k)').^2); % in m/s
% erreur_vz = sqrt((vzm(1:k) - stateEst(6,1:k)').^2); % in m/s
% 
% % error for the state without corrections
% erreur_xNoCorr = sqrt((xm(1:k) - stateNoCorr(1,1:k)').^2); % in m
% erreur_yNoCorr = sqrt((ym(1:k) - stateNoCorr(2,1:k)').^2); % in m
% erreur_zNoCorr = sqrt((zm(1:k) - stateNoCorr(3,1:k)').^2); % in m
% erreur_vxNoCorr = sqrt((vxm(1:k) - stateNoCorr(4,1:k)').^2); % in m/s
% erreur_vyNoCorr = sqrt((vym(1:k) - stateNoCorr(5,1:k)').^2); % in m/s
% erreur_vzNoCorr = sqrt((vzm(1:k) - stateNoCorr(6,1:k)').^2); % in m/s
% 
% % error for the state no ekf
% erreur_xNoEkf = sqrt((xm(1:k) - stateNoEkf(1,1:k)').^2); % in m
% erreur_yNoEkf  = sqrt((ym(1:k) - stateNoEkf(2,1:k)').^2); % in m
% erreur_zNoEkf  = sqrt((zm(1:k) - stateNoEkf(3,1:k)').^2); % in m
% erreur_vxNoEkf  = sqrt((vxm(1:k) - stateNoEkf(4,1:k)').^2); % in m/s
% erreur_vyNoEkf  = sqrt((vym(1:k) - stateNoEkf(5,1:k)').^2); % in m/s
% erreur_vzNoEkf  = sqrt((vzm(1:k) - stateNoEkf(6,1:k)').^2); % in m/s
% 
% % plots
% xGraph  = 1:1:k;
% 
% figure
% subplot(2,2,1)
% plot(xGraph,xm(1:k),xGraph,xm(1:k),'ob',xGraph,stateEst(1,1:k),'*r')
% title('Subplot 1: xm')
% subplot(2,2,2)
% plot(xGraph,ym(1:k),xGraph,ym(1:k),'ob',xGraph,stateEst(2,1:k),'*r')
% title('Subplot 2: ym')
% subplot(2,2,3)
% plot(xGraph,zm(1:k),xGraph,zm(1:k),'ob',xGraph,stateEst(3,1:k),'*r')
% title('Subplot 3: zm')
% 
% 
% X_desired = 0:time_vector(end);
% 
% figure
% subplot(2,2,1)
% Y1_desired = interp1(time_vector, xm(1:k), X_desired);
% Y2_desired = interp1(time_vector, stateEst(1,1:k), X_desired);
% Y3_desired = interp1(time_vector, stateNoCorr(1,1:k), X_desired);
% plot(X_desired,Y1_desired,X_desired,Y2_desired,X_desired,Y3_desired)
% xlabel('time (s)')
% ylabel('position (m)')
% title('x')
% subplot(2,2,2)
% Y1_desired = interp1(time_vector, ym(1:k), X_desired);
% Y2_desired = interp1(time_vector, stateEst(2,1:k), X_desired);
% Y3_desired = interp1(time_vector, stateNoCorr(2,1:k), X_desired);
% plot(X_desired,Y1_desired,X_desired,Y2_desired,X_desired,Y3_desired)
% xlabel('time (s)')
% ylabel('position (m)')
% title('y')
% subplot(2,2,3)
% Y1_desired = interp1(time_vector, zm(1:k), X_desired);
% Y2_desired = interp1(time_vector, stateEst(3,1:k), X_desired);
% Y3_desired = interp1(time_vector, stateNoCorr(3,1:k), X_desired);
% plot(X_desired,Y1_desired,X_desired,Y2_desired,X_desired,Y3_desired)
% xlabel('time (s)')
% ylabel('position (m)')
% title('z')
% legend("GPS","Estimate with corrections","Estimate wihtout corrections");
% 
% 
% figure
% subplot(2,2,1)
% Y1_desired = interp1(time_vector, vxm(1:k), X_desired);
% Y2_desired = interp1(time_vector, stateEst(4,1:k), X_desired);
% Y3_desired = interp1(time_vector, stateNoCorr(4,1:k), X_desired);
% plot(X_desired,Y1_desired,X_desired,Y2_desired,X_desired,Y3_desired)
% xlabel('time (s)')
% ylabel('velocity (m/s)')
% title('vx')
% subplot(2,2,2)
% Y1_desired = interp1(time_vector, vym(1:k), X_desired);
% Y2_desired = interp1(time_vector, stateEst(5,1:k), X_desired);
% Y3_desired = interp1(time_vector, stateNoCorr(5,1:k), X_desired);
% plot(X_desired,Y1_desired,X_desired,Y2_desired,X_desired,Y3_desired)
% xlabel('time (s)')
% ylabel('velocity (m/s)')
% title('vy')
% subplot(2,2,3)
% Y1_desired = interp1(time_vector, vzm(1:k), X_desired);
% Y2_desired = interp1(time_vector, stateEst(6,1:k), X_desired);
% Y3_desired = interp1(time_vector, stateNoCorr(6,1:k), X_desired);
% plot(X_desired,Y1_desired,X_desired,Y2_desired,X_desired,Y3_desired)
% xlabel('time (s)')
% ylabel('velocity (m/s)')
% title('vz')
% legend("GPS","Estimate with corrections","Estimate wihtout corrections");
% 
% 
% 
% figure
% subplot(2,2,1)
% Y1_desired = interp1(time_vector, xm(1:k), X_desired);
% Y2_desired = interp1(time_vector, stateEst(1,1:k), X_desired);
% Y3_desired = interp1(time_vector, stateNoCorr(1,1:k), X_desired);
% Y4_desired = interp1(time_vector, stateNoEkf(1,1:k), X_desired);
% plot(X_desired,Y1_desired,X_desired,Y2_desired,X_desired,Y3_desired,X_desired,Y4_desired)
% xlabel('time (s)')
% ylabel('position (m)')
% title('x')
% subplot(2,2,2)
% Y1_desired = interp1(time_vector, ym(1:k), X_desired);
% Y2_desired = interp1(time_vector, stateEst(2,1:k), X_desired);
% Y3_desired = interp1(time_vector, stateNoCorr(2,1:k), X_desired);
% Y4_desired = interp1(time_vector, stateNoEkf(2,1:k), X_desired);
% plot(X_desired,Y1_desired,X_desired,Y2_desired,X_desired,Y3_desired,X_desired,Y4_desired)
% xlabel('time (s)')
% ylabel('position (m)')
% title('y')
% subplot(2,2,3)
% Y1_desired = interp1(time_vector, zm(1:k), X_desired);
% Y2_desired = interp1(time_vector, stateEst(3,1:k), X_desired);
% Y3_desired = interp1(time_vector, stateNoCorr(3,1:k), X_desired);
% Y4_desired = interp1(time_vector, stateNoEkf(3,1:k), X_desired);
% plot(X_desired,Y1_desired,X_desired,Y2_desired,X_desired,Y3_desired,X_desired,Y4_desired)
% xlabel('time (s)')
% ylabel('position (m)')
% title('z')
% legend("GPS","Estimate with corrections","Estimate wihtout corrections","nothing");