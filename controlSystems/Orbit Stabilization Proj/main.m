% the purpose of this script is to setup the framework to complete the
% non-linear systems project. We're looking here to design a controller to
% stabilise the orbit to a given set of parameters. Here is where all the
% heavy lifting for the project is done. As well, it will display all of
% the graphs and visuals to move forwards. Recall, we have two distince
% reference frames here, which means the angular velocity will be in one of
% them and the linear position in the other. 
close all
clear all


% we need here the moments of inertial for this spacecraft. We're going to
% assume that this is much simpler than it actually is for the scope of
% this project, nice and easy shapes for moment of inertia.   
spacecraft_MoI_N = 0.17; % kg m^2, about the 10 cm axis
spacecraft_MoI_E = 0.14; % kg m^2, about the 20 cm axis
spacecraft_MoI_D = 0.05; % kg m^2, about the 36 cm axis
reactionWheel_MoI = 0.0023; % kg m^2

% some physical parameters as reference data to run the simulation. 
mu = 3.986004418e14; % standart gravitational parameter (m^3/s^2)
R = 6.37813e6; % Earth radius (m)
J2 = 1.082626e-3; % dimensional coefficient

% angular velocity of the reference frame. We need this guy and it can be
% quickly obtained from a calculation 


timeSpan = [0 86400]; % seconds in a day 
initialPosition = [0.1,pi/4,pi/5];
[ww,xx,yy,zz] = euler2quat(initialPosition(1),initialPosition(2),initialPosition(3),'zyx');
quatVect = [ww;xx;yy;zz];
y0 = [-549214,2450301,6695830,4463,-5485,2383,ww,xx,yy,zz,0,0,0,0,0,0]; % initial state vector
refFrameAngVel = cross(y0(1:3),y0(4:6))./(norm(y0(1:3)).^2);  
y0(11:13) = refFrameAngVel;
referenceTrajectory = []; % set of orbital parameters used stabilize the orbit. 

% this is the line of code to run the ode solver. We have multiple choices
% of ode solvers in matlab. Look through documentation and choose
% accordingly. 
[t,y] = ode78(@(t,y) findStateDot(t,y,referenceTrajectory,spacecraft_MoI_N,...
    spacecraft_MoI_E, spacecraft_MoI_D, reactionWheel_MoI, mu, R,J2),timeSpan,y0);


% we now want to plot the outputs here
figure
subplot(2,3,1)
plot(t,y(:,1))
title('X Position as a Function of Time')
subplot(2,3,2)
plot(t,y(:,2))
title('Y Position as a Function of Time')
subplot(2,3,3)
plot(t,y(:,3))
title('Z Position as a Function of Time')
subplot(2,3,4)
plot(t,y(:,4))
title('X Velocity as a Function of Time')
subplot(2,3,5)
plot(t,y(:,5))
title('Y Velocity as a Function of Time')
subplot(2,3,6)
plot(t,y(:,6))
title('Z Velocity as a Function of Time')
sgtitle('Linear Position and Velocity for the Spacecraft')



for i = 1:length(y)
    [roll(i), pitch(i), yaw(i)] = q2e(y(i,7),y(i,8),y(i,9),y(i,10));
end
figure
plot(t,[roll',pitch',yaw'])
legend('roll','pitch','yaw')

% we now want the angular output plots here
figure
subplot(2,3,1)
plot(t,y(:,14))
title('Reaction Wheel Angular Velocity about N as a Function of Time')
subplot(2,3,2)
plot(t,y(:,15))
title('Reaction Wheel Angular Velocity about E as a Function of Time')
subplot(2,3,3)
plot(t,y(:,16))
title('Reaction Wheel Angular Velocity about D as a Function of Time')
subplot(2,3,4)
plot(t,y(:,11))
title('Angular Velocity about N as a Function of Time')
subplot(2,3,5)
plot(t,y(:,12))
title('Angular Velocity about E as a Function of Time')
subplot(2,3,6)
plot(t,y(:,13))
title('Angular Velcoity about D as a Function of Time')
sgtitle('Angular Velocities for Spacecraft and Reaction Wheels')






