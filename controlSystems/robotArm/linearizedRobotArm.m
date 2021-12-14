function xdot = linearizedRobotArm(t,y,A,B,C,controlInp)
% the purpose of this function is to map out the linearized version of the
% robotic arm. We're going to compare how the linearized version does under
% perturbation compared to the first one.
xdot = (A * y) + (B .* controlInp);
