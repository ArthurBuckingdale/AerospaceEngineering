function controlInput = computeControlInput()
% the purpose of this function is to compute the control input to stabilize
% the spacecraft. It will take in the current state, then the reference 
% trajectory to compute the control vector u which will enter into the
% propagator. 