function dot = findStateDot(t,state,referenceTrajectory,spacecraft_MoI_N,...
    spacecraft_MoI_E, spacecraft_MoI_D, reactionWheel_MoI, mu, R,J2)
% the purpose of this function is to propagate the spacecraft. It will take
% in the state, and apply the force model to the linear and angular
% positions of the state. What's very important to note here, is how the
% angular velocity of the spacecraft and the reaction wheels are present
% here. The reaction wheel angular velocity affects the overall spacecraft
% angular velocity; Also note, that we needed to attach the controller
% variables to this. 

% first up, we have the state variables being input here. 
x = state(1);   %(m)
y = state(2);   %(m)
z = state(3);   %(m)
vx = state(4);  %(m/s)
vy = state(5);  %(m/s)
vz = state(6);  %(m/s)
phiN = state(14); % rad/s reaction wheel angular speed
phiE = state(15); % rad/s reaction wheel angular speed
phiD = state(16); % rad/s reaction wheel angular speed
omegaN = state(11) - (phiN * (reactionWheel_MoI/spacecraft_MoI_N)); %rad/s angular speed
omegaE = state(12) - (phiE * (reactionWheel_MoI/spacecraft_MoI_E)); %rad/s angular speed
omegaD = state(13) - (phiD * (reactionWheel_MoI/spacecraft_MoI_D)); %rad/s angular speed


% this is where we use the reference trajectory to compute the necesary
% control input. To simplify this problem, I am going to make the
% trajectory a vector of something like [altitude, velocity, ang Position].
% I don't want to program in a whole set of state vectors as reference
% trajectories here. That just seems like a huge pain in the ass. Even
% though I have state vectors from NEOSSat.

% controlInput = computeControlInput(state, referenceTrajectory);
controlInput = zeros(7);
% now, we have the control variables from the input, it is labelled as
% state, but that is simply just to pass it in here.
ionThrust = controlInput(1);%controlInput(1); % linear ion thruster
tRN = 0.0000001;%controlInput(2);  % torque of reaction wheel about N axis
tMtN = controlInput(3); % torque of magnetorque rod about N axis
tRE = controlInput(4);  % torque of reaction wheel about E axis
tMtE = controlInput(5); % torque of magnetorque rod about E axis
tRD = controlInput(6);  % torque of reaction wheel about D axis
tMtD = controlInput(7); % torque of magnetorque rod about D axis

% first step, get out angular velocities in the state they need to be in
angularVelocityMatrix(2:4,1) = [omegaN,omegaE,omegaD];
angularVelocityMatrix(1,2:4) = -[omegaN,omegaE,omegaD];
angularVelocityMatrix(2,2:end) = cross([omegaN,omegaE,omegaD]',[1,0,0]');
angularVelocityMatrix(3,2:end) = cross([omegaN,omegaE,omegaD]',[0,1,0]');
angularVelocityMatrix(4,2:end) = cross([omegaN,omegaE,omegaD]',[0,0,1]');

anglesDeltas = (((1/2) * angularVelocityMatrix * state(7:10)));

% 
% [rollDelta, pitchDelta, yawDelta] = q2e(anglesDeltas(1),anglesDeltas(2),anglesDeltas(3),anglesDeltas(4));

r = norm([x,y,z]); % norm between center of Earth and the satellite position

dot = zeros(15,1);

dot(1) = vx;
dot(2) = vy;
dot(3) = vz;
dot(4) = -(mu*x)/r^3 + ((mu*J2*R^2)/2)*(((15*x*z^2)/r^7) - (3*x)/r^5) + ionThrust;
dot(5) = -(mu*y)/r^3 + ((mu*J2*R^2)/2)*(((15*y*z^2)/r^7) - (3*y)/r^5);
dot(6) = -(mu*z)/r^3 + ((mu*J2*R^2)/2)*(((15*z^3)/r^7) - (9*z)/r^5);
dot(7) = anglesDeltas(1); % I need to add the angular momentum difference here
dot(8) = anglesDeltas(2);
dot(9) = anglesDeltas(3);
dot(10) = anglesDeltas(4);
dot(11) = (tRN + tMtN) / (spacecraft_MoI_N + reactionWheel_MoI);
dot(12) = (tRE + tMtE) / (spacecraft_MoI_E + reactionWheel_MoI);
dot(13) = (tRD + tMtD) / (spacecraft_MoI_D + reactionWheel_MoI);
dot(14) = tRN / (reactionWheel_MoI);
dot(15) = tRE / (reactionWheel_MoI);
dot(16) = tRD / (reactionWheel_MoI);
end











