function xDot = linearizedBioReactor(t,x,A)
% the purpose of this function is to give the linearized bio reactor.
xDot = A * x;