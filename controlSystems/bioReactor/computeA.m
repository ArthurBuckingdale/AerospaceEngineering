function A = computeA(x,D)
% the purpose of this function is to compute the linearized A matrix at a
% given point and then use it for the linearize function.
muBar = 0.53;
kM = 0.12;
k1 = 0.5;
Z = 0.5;

A(1,1) = muBioReactor(x,muBar,kM,k1) - D;
A(1,2) = ((muBioReactor(x,muBar,kM,k1) - muBioReactor(x+0.0001,muBar,kM,k1))/0.0001) * x(1);
A(2,1) = -1/Z * muBioReactor(x,muBar,kM,k1);
A(2,2) = -D - (1/Z) * x(1) * ((muBioReactor(x,muBar,kM,k1) - muBioReactor(x+0.0001,muBar,kM,k1))/0.0001);