function xDot = bioReactor(t,x,D,X2f)
% the purpose of this function is to simulate the differential equations of
% "motion" for the bio reactor. Pretty cool

%declare some constants given by the question. Presumably these could
%change, depending on the system I presume. 
muBar = 0.53;
kM = 0.12;
k1 = 0.5;
Z = 0.5;

xDot(1,1) = (muBioReactor(x,muBar,kM,k1) - D) * x(1);
xDot(2,1) = D * (X2f - x(2)) - (1/Z) * x(1) * muBioReactor(x,muBar,kM,k1);
