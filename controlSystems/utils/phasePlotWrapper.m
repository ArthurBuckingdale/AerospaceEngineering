% the purpose of this script is to handle the filling of the phase plot
% function.
close all
clear all

% we want to create a list of initial conditions to draw the phase plot
xlimit = [-10,10];
ylimit = [-10,10];
n=1;
for i = xlimit(1):2:xlimit(2)
    for j = ylimit(1):2:ylimit(2)
        initialCond(n,:) = [i,j];
        n=n+1;
    end
end

timeSpan = [0,5];
funHandle = @MID1Q2;
phasePlot(initialCond, timeSpan, funHandle)