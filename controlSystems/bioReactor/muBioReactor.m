function output = muBioReactor(x,muBar,kM,k1)
% the purpose of this function is to provide a calculation for the embedded
% function in the bioreactor question.
output = (muBar*x(2))./(kM + x(2) + k1*x(2).^2);
