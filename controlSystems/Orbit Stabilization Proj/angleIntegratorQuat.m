function [roll, pitch, yaw] = angleIntegratorQuat(angularState, timeDelta)
% the purpose of this function is to integrate the quaternion angles. using
% the formulae from navigaation systems homework.

[w,x,y,z] = euler2quat(angularState(1),angularState(2),angularState(3),'zyx');
quatVect = [w;x;y;z];

% first step, get out angular velocities in the state they need to be in
angularVelocityMatrix(2:4,1) = angularState(4:6);
angularVelocityMatrix(1,2:4) = -angularState(4:6);
angularVelocityMatrix(2,2:end) = cross(angularState(4:6)',[1,0,0]');
angularVelocityMatrix(3,2:end) = cross(angularState(4:6)',[0,1,0]');
angularVelocityMatrix(4,2:end) = cross(angularState(4:6)',[0,0,1]');

angles = ((cos((timeDelta/2) * norm(angularState(4:6))) * eye(4)) + ...
    sin((timeDelta/2) * norm(angularState(4:6)))/norm(angularState(4:6)+0.0000000000000001) * angularVelocityMatrix) * quatVect;

[roll, pitch, yaw] = q2e(angles(1),angles(2),angles(3),angles(4));