function position = robotArm(t,y,controlInp)
% the purpose of this function is to simulate the robotic arm from a
% homework question.
position = [y(2); -10 * sin(y(1)) - y(2) + controlInp]