clc
clear
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
sigma_velocity = 0.1; %(m/s)
sigma_measure = 50; % (m/s)
sigma_position = 100; %meters

%generate noise
W = sigma_velocity*randn(1,length(xm));
V = sigma_measure*randn(1,length(xm));

%matrices
G = [0 0 0;
    0 0 0;
    0 0 0;
    1 0 0;
    0 1 0;
    0 0 1;];

D = [1 0 0;
    0 1 0;
    0 0 1];

P = [sigma_position^2 0 0 0 0 0;
    0 sigma_position^2 0 0 0 0;
    0 0 sigma_position^2 0 0 0;
    0 0 0 sigma_velocity^2 0 0;
    0 0 0 0 sigma_velocity^2 0;
    0 0 0 0 0 sigma_velocity^2;]; % initial guess of P

C = [0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1;];

stateEst = zeros(6,length(xm)); % length will probably change when we merge code
x0 = [xm(1);ym(1);zm(1);vxm(1);vym(1);vzm(1)]; % initiale value
stateEst(:,1) = x0;
stateNoCorr = stateEst; % the state without any corrections
stateNoEkf = stateEst;

trueState = [xm,ym,zm,vxm,vym,vzm]; % true measurements

interval = zeros(length(xm),1); % time interval between measurements


% EKF
Kmax = 24; % %%%%%%%%%%%_MODIFY__%%%%%%%%%%%% change this if you want to check more or less values, Max is 100;

time_vector = zeros(Kmax,1);

for k = 2:Kmax
    
    interval(k) = time(k) - time(k-1);
    time_vector(k) = time(k) - time(1);
    
    % prediction
    [~,integratedState] = ode45(@findStateDot, [0 interval(k)], stateEst(:,k-1));
    X_moins = integratedState(end,:)'; % state prediction
    stateNoCorr(:,k) = X_moins; % for comparaison only, not neccessary
    F = findJacobian(X_moins);
    [~,integratedCov] = ode45(@(~,integratedCov) findCovarianceDot(P,F,G,W(k)), [0 interval(k)],P);
    P_moins = reshape(integratedCov(end,:),6,6); % covariance predecition
    
    % measurement update
    K = P_moins*C'/(C*P_moins*C' + D*V(k)*D');
    velocity_noise = 50*randn(3,1);
    X_plus = X_moins + K*(trueState(k,4:6)' + velocity_noise - C*X_moins); % update state
    stateEst(:,k) = X_plus;
    P = (eye(6) - K*C)*P_moins; % update covariance estimate
end


for k = 2:Kmax % we correct the value only when new measurement is available
    
    [~,integratedState] = ode45(@findStateDot, [0 interval(k)], stateEst(:,k-1));
    X_moins = integratedState(end,:)';
    
    stateNoEkf(:,k) = X_moins;
    
end


% error for the corrected state
erreur_x = sqrt((xm(1:k) - stateEst(1,1:k)').^2); % in m
erreur_y = sqrt((ym(1:k) - stateEst(2,1:k)').^2); % in m
erreur_z = sqrt((zm(1:k) - stateEst(3,1:k)').^2); % in m
erreur_vx = sqrt((vxm(1:k) - stateEst(4,1:k)').^2); % in m/s
erreur_vy = sqrt((vym(1:k) - stateEst(5,1:k)').^2); % in m/s
erreur_vz = sqrt((vzm(1:k) - stateEst(6,1:k)').^2); % in m/s

% error for the state without corrections
erreur_xNoCorr = sqrt((xm(1:k) - stateNoCorr(1,1:k)').^2); % in m
erreur_yNoCorr = sqrt((ym(1:k) - stateNoCorr(2,1:k)').^2); % in m
erreur_zNoCorr = sqrt((zm(1:k) - stateNoCorr(3,1:k)').^2); % in m
erreur_vxNoCorr = sqrt((vxm(1:k) - stateNoCorr(4,1:k)').^2); % in m/s
erreur_vyNoCorr = sqrt((vym(1:k) - stateNoCorr(5,1:k)').^2); % in m/s
erreur_vzNoCorr = sqrt((vzm(1:k) - stateNoCorr(6,1:k)').^2); % in m/s

% error for the state no ekf
erreur_xNoEkf = sqrt((xm(1:k) - stateNoEkf(1,1:k)').^2); % in m
erreur_yNoEkf  = sqrt((ym(1:k) - stateNoEkf(2,1:k)').^2); % in m
erreur_zNoEkf  = sqrt((zm(1:k) - stateNoEkf(3,1:k)').^2); % in m
erreur_vxNoEkf  = sqrt((vxm(1:k) - stateNoEkf(4,1:k)').^2); % in m/s
erreur_vyNoEkf  = sqrt((vym(1:k) - stateNoEkf(5,1:k)').^2); % in m/s
erreur_vzNoEkf  = sqrt((vzm(1:k) - stateNoEkf(6,1:k)').^2); % in m/s

% plots
xGraph  = 1:1:k;

figure
subplot(1,3,1)
plot(xGraph,xm(1:k),'-ob',xGraph,stateEst(1,1:k),'*r')
title('GPS versus Propagated Coordinates in X dim')
xlabel('Time (hr)')
ylabel('Postion in X (m)')
legend('GPS','Propagated State')
subplot(1,3,2)
plot(xGraph,ym(1:k),'-ob',xGraph,stateEst(2,1:k),'*r')
title('GPS versus Propagated Coordinates in Y dim')
xlabel('Time (hr)')
ylabel('Postion in Y (m)')
legend('GPS','Propagated State')
subplot(1,3,3)
plot(xGraph,zm(1:k),'-ob',xGraph,stateEst(3,1:k),'*r')
title('GPS versus Propagated Coordinates in Z dim')
xlabel('Time (hr)')
ylabel('Postion in Z (m)')
legend('GPS','Propagated State')
sgtitle('Propagated versus GPS Coordinates for NEOSSat')

X_desired = 0:time_vector(end);

figure
subplot(2,2,1)
Y1_desired = interp1(time_vector, xm(1:k), X_desired);
Y2_desired = interp1(time_vector, stateEst(1,1:k), X_desired);
Y3_desired = interp1(time_vector, stateNoCorr(1,1:k), X_desired);
plot(X_desired,Y1_desired,X_desired,Y2_desired,X_desired,Y3_desired)
xlabel('time (s)')
ylabel('position (m)')
title('x')
subplot(2,2,2)
Y1_desired = interp1(time_vector, ym(1:k), X_desired);
Y2_desired = interp1(time_vector, stateEst(2,1:k), X_desired);
Y3_desired = interp1(time_vector, stateNoCorr(2,1:k), X_desired);
plot(X_desired,Y1_desired,X_desired,Y2_desired,X_desired,Y3_desired)
xlabel('time (s)')
ylabel('position (m)')
title('y')
subplot(2,2,3)
Y1_desired = interp1(time_vector, zm(1:k), X_desired);
Y2_desired = interp1(time_vector, stateEst(3,1:k), X_desired);
Y3_desired = interp1(time_vector, stateNoCorr(3,1:k), X_desired);
plot(X_desired,Y1_desired,X_desired,Y2_desired,X_desired,Y3_desired)
xlabel('time (s)')
ylabel('position (m)')
title('z')
legend("GPS","Estimate with corrections","Estimate wihtout corrections");


figure
subplot(2,2,1)
Y1_desired = interp1(time_vector, vxm(1:k), X_desired);
Y2_desired = interp1(time_vector, stateEst(4,1:k), X_desired);
Y3_desired = interp1(time_vector, stateNoCorr(4,1:k), X_desired);
plot(X_desired,Y1_desired,X_desired,Y2_desired,X_desired,Y3_desired)
xlabel('time (s)')
ylabel('velocity (m/s)')
title('vx')
subplot(2,2,2)
Y1_desired = interp1(time_vector, vym(1:k), X_desired);
Y2_desired = interp1(time_vector, stateEst(5,1:k), X_desired);
Y3_desired = interp1(time_vector, stateNoCorr(5,1:k), X_desired);
plot(X_desired,Y1_desired,X_desired,Y2_desired,X_desired,Y3_desired)
xlabel('time (s)')
ylabel('velocity (m/s)')
title('vy')
subplot(2,2,3)
Y1_desired = interp1(time_vector, vzm(1:k), X_desired);
Y2_desired = interp1(time_vector, stateEst(6,1:k), X_desired);
Y3_desired = interp1(time_vector, stateNoCorr(6,1:k), X_desired);
plot(X_desired,Y1_desired,X_desired,Y2_desired,X_desired,Y3_desired)
xlabel('time (s)')
ylabel('velocity (m/s)')
title('vz')
legend("GPS","Estimate with corrections","Estimate wihtout corrections");



figure
subplot(2,2,1)
Y1_desired = interp1(time_vector, xm(1:k), X_desired);
Y2_desired = interp1(time_vector, stateEst(1,1:k), X_desired);
Y3_desired = interp1(time_vector, stateNoCorr(1,1:k), X_desired);
Y4_desired = interp1(time_vector, stateNoEkf(1,1:k), X_desired);
plot(X_desired,Y1_desired,X_desired,Y2_desired,X_desired,Y3_desired,X_desired,Y4_desired)
xlabel('time (s)')
ylabel('position (m)')
title('x')
subplot(2,2,2)
Y1_desired = interp1(time_vector, ym(1:k), X_desired);
Y2_desired = interp1(time_vector, stateEst(2,1:k), X_desired);
Y3_desired = interp1(time_vector, stateNoCorr(2,1:k), X_desired);
Y4_desired = interp1(time_vector, stateNoEkf(2,1:k), X_desired);
plot(X_desired,Y1_desired,X_desired,Y2_desired,X_desired,Y3_desired,X_desired,Y4_desired)
xlabel('time (s)')
ylabel('position (m)')
title('y')
subplot(2,2,3)
Y1_desired = interp1(time_vector, zm(1:k), X_desired);
Y2_desired = interp1(time_vector, stateEst(3,1:k), X_desired);
Y3_desired = interp1(time_vector, stateNoCorr(3,1:k), X_desired);
Y4_desired = interp1(time_vector, stateNoEkf(3,1:k), X_desired);
plot(X_desired,Y1_desired,X_desired,Y2_desired,X_desired,Y3_desired,X_desired,Y4_desired)
xlabel('time (s)')
ylabel('position (m)')
title('z')
legend("GPS","Estimate with corrections","Estimate wihtout corrections","nothing");