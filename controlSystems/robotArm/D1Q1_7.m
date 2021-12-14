% the purpose of this script is to answer part D of the question 1.7 from
% the systemes multivariables homework.
close all
clear all

% time duration
timeSpan = [0 10];

% initial position of the robotic arm
y0 = [3 * pi/4; 0;];

% control input
controlInp = 5*sqrt(2);


% this is the section where we are generating the ground truth using the
% matlab ode function. We will then perturb our linearize model and
% propagate it to see what happens as the perturbation gets larger.
% running the ode
[t,y] = ode45(@(t,y) robotArm(t,y,controlInp),timeSpan,y0);

% plotting of this function
figure
plot(t,y,'*')
title('Angular Position and Velocity as a Function of Time Ground Truth')
ylabel('Angular Position and Velocity')
xlabel('Time (s)')
legend('Pos','Vel')


% we now want to plug the the linearized robot arm function into a loop and
% see what happens as we slightly perturb the control input
A = [0,1;5*sqrt(2),-1];
B = [0;1];
C = [1, 0];
thetaEquil = 3 * pi/4;
controlInp = 0.001;
x(1,:) = [thetaEquil,0];
t=0;
x0(:,1) = [0;0];
%compute first order linear correction

figure
hold on
for j = -0.01:0.002:0.01
    for i = 1:50
        linearizedRobotArm(t,[0;0],A,B,C,controlInp)
        deltas = linearizedRobotArm(t,x0(:,i),A,B,C,j) * 0.1;
        x0(:,i+1) = x0(:,i) + deltas;
        x(i+1,:) = x(i,:) + deltas';
    end
    plot(x(:,1),'*')
end
hold off
title('Linearized Propagated System')
xlabel('Time')
ylabel('Position')



% let's see what happens when we perturb the non linear formula
m=1;
figure
hold on
y0 = [3 * pi/4; 0;];
controlInp = 5*sqrt(2);
for i = -0.04:0.01:0.04
    [t3(m).val,y3(m).val] = ode45(@(t,y) robotArm(t,y,controlInp + i),timeSpan,y0);
    plot(t3(m).val,y3(m).val(:,1),'*')
    m=m+1;
end
hold off






