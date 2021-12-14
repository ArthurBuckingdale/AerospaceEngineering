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
time = posixtime(D{:,9});

stateEst = zeros(6,length(xm)); % length will probably change when we merge code
x0 = [xm(1);ym(1);zm(1);vxm(1);vym(1);vzm(1)]; % initiale value
stateEst(:,1) = x0;

trueState = [xm,ym,zm,vxm,vym,vzm]; % true measurements

interval = zeros(length(xm),1); % time interval between tru measurements

% begin loop
for k = 2:length(xm) % we correct the value only when new measurement is available
    
    interval(k) = time(k) - time(k-1); % find the time interval between measurements
   
    [t,y] = ode45(@findStateDot, [0 interval(k)], trueState(k-1,:));
    stateEst(:,k) = y(end,:)';
end

% computes energy
mu = 3.986004418e14; % standart gravitational parameter (m^3/s^2)
E_est = zeros(k,1);
E_reel = zeros(k,1);

for i = 1:k
    E_est(i) = norm([stateEst(4,i),stateEst(5,i),stateEst(6,i)])^2/2 - mu/norm([stateEst(1,i),stateEst(2,i),stateEst(3,i)]);
    E_reel(i) = norm([vxm(i),vym(i),vzm(i)])^2/2 - mu/norm([xm(i),ym(i),zm(i)]);
end


% computes error 
erreur_x = sqrt((xm(1:k) - stateEst(1,1:k)').^2); % in m
erreur_y = sqrt((ym(1:k) - stateEst(2,1:k)').^2); % in m
erreur_z = sqrt((zm(1:k) - stateEst(3,1:k)').^2); % in m
erreur_energy = sqrt((E_est - E_reel).^2); % m^2/s^2

% plots
xGraph  = 1:1:length(xm);

figure
subplot(2,2,1)
plot(xGraph,xm,xGraph,xm,'ob',xGraph,stateEst(1,1:k),'*r')
xlabel('sample number') 
ylabel('position (m)') 
title('x')
subplot(2,2,2)
plot(xGraph,ym,xGraph,ym,'ob',xGraph,stateEst(2,1:k),'*r')
xlabel('sample number') 
ylabel('position (m)') 
title('y')
subplot(2,2,3)
plot(xGraph,zm,xGraph,zm,'ob',xGraph,stateEst(3,1:k),'*r')
xlabel('sample number') 
ylabel('position (m)') 
title('z')


figure
subplot(2,2,1)
plot(xGraph,erreur_x,'ob',xGraph,erreur_energy,'*r')
xlabel('sample number') 
ylabel('energy') 
title('error in x')
subplot(2,2,2)
plot(xGraph,erreur_y,'ob',xGraph,erreur_energy,'*r')
xlabel('sample number') 
ylabel('energy') 
title('error in y')
subplot(2,2,3)
plot(xGraph,erreur_z,'ob',xGraph,erreur_energy,'*r')
xlabel('sample number') 
ylabel('energy') 
title('error in z')
legend("position error","energy error");