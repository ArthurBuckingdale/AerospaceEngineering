% the purpose of this script is to solve question 2 part three part of
% assignment three for the navigation systems class.
close all
clear all

% we first off need to define an angular trajectory for each different of
% the attitude angles. Our tests require 720 degrees rotation about each
% axis.
vect = -180:1:180;
vect2 = -44:1:44;

psiInitial =   0.1*ones(1,728);%vect3;%
%psiInitial = psiInitial + 0.05*randn(1,length(psiInitial));
thetaInitial = [sind(vect2),sind(fliplr(vect2)),sind(vect2),sind(fliplr(vect2))];%zeros(1,722);%
thetaInitial = 0.5*[thetaInitial,thetaInitial];% + 0.01*randn(1,2*length(thetaInitial));
phiInitial = sind([vect,vect]); %zeros(1,728);%
%phiInitial = phiInitial + 0.05*randn(1,length(phiInitial));

figure
subplot(1,3,1)
plot(psiInitial)
title('\psi Trajectory')
xlabel('Time')
ylabel('Angle')
subplot(1,3,2)
plot(thetaInitial)
title('\theta Trajectory')
xlabel('Time')
ylabel('Angle')
subplot(1,3,3)
plot(phiInitial)
title('\phi Trajectory')
xlabel('Time')
ylabel('Angle')
hold off


% we now need to convert these back into a rotation matrix.
for i = 1:length(thetaInitial)
    rotationMatrix(1,1,i) = cosd(thetaInitial(i)) .* cosd(psiInitial(i));
    rotationMatrix(1,2,i) = -cosd(phiInitial(i)).*sind(psiInitial(i)) +...
        sind(phiInitial(i)).*sind(thetaInitial(i)).*cosd(psiInitial(i));
    rotationMatrix(1,3,i) = sind(phiInitial(i)).*sind(psiInitial(i)) + ...
        cosd(phiInitial(i)).*sind(thetaInitial(i)).*cosd(psiInitial(i));
    rotationMatrix(2,1,i) = cosd(thetaInitial(i)).*sind(psiInitial(i));
    rotationMatrix(2,2,i) = cosd(phiInitial(i)).*cosd(psiInitial(i)) + ...
        sind(phiInitial(i)).*sind(thetaInitial(i)).*sind(psiInitial(i));
    rotationMatrix(2,3,i) = -sind(phiInitial(i)).*cos(psiInitial(i)) + ...
        cosd(phiInitial(i)).*sind(thetaInitial(i)).*sind(psiInitial(i));
    rotationMatrix(3,1,i) = -sind(thetaInitial(i));
    rotationMatrix(3,2,i) = sind(phiInitial(i)).*cosd(thetaInitial(i));
    rotationMatrix(3,3,i) = cosd(phiInitial(i)).*cosd(thetaInitial(i));
end

% we must now numerically differentiate the rotation matrix
for i = 1:size(rotationMatrix,3)-1
    %assuming 1 deg/sec
    rotationMatrixDot(:,:,i) = rotationMatrix(:,:,i+1) - rotationMatrix(:,:,i);
end

% we now determine the \Omega matrix for each and pick out each component
for i = 1:size(rotationMatrix,3)-1
    omegaMatrix(:,:,i) = rotationMatrix(:,:,i) \ rotationMatrixDot(:,:,i);
    omegaMatrix(:,:,i) = omegaMatrix(:,:,i) - diag(diag(omegaMatrix(:,:,i)));
    angularVelocity(i,1) = omegaMatrix(2,3,i) + 0.00005*randn(1,1);
    angularVelocity(i,2) = omegaMatrix(3,1,i) + 0.00005*randn(1,1);
    angularVelocity(i,3) = omegaMatrix(1,2,i) + 0.00005*randn(1,1);
end

% last but not least, let's plot these recovered angular velocities
figure
hold on
plot(angularVelocity(:,1))
plot(angularVelocity(:,2))
plot(angularVelocity(:,3))
hold off
title('Computed Angular Velocities')
xlabel('Time')
legend('\omega_1','\omega_2','\omega_3')

% this is the money ticket where we see how good our integrator is.
initialRotationMatrix = rotationMatrix(:,:,1);
timeStep = 1;
roll(1) = pi/3;
pitch(1) = pi/4;
yaw(1) = 0.1;
numObs = 40000;
for i = 1:numObs
    [roll(i+1),pitch(i+1),yaw(i+1)] = angleIntegratorQuat([roll(i), pitch(i), yaw(i),0.01,0.01,0.01], 0.01);
end


timeStamps = 0:0.01:(0.01 * numObs);

figure
subplot(2,2,1)
plot(timeStamps, roll,'-*')
subplot(2,2,2)
plot(timeStamps, pitch,'-*')
subplot(2,2,3)
plot(timeStamps, yaw,'-*')

























