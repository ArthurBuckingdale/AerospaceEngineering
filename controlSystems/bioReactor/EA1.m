% the purpose of this script is to answer the questions from the additional
% exercises in multivariable systems. The system is a bioreactor or
% fermentation device to make beer for example. 
clear all
close all

%% using the trim algo for part A
% we're looking for the roots of a non-linear system to see if i correctly
% calculated my equilibria. Can use the trim furnction but I don't have it
% because matlab.
D = 0.3;
X2f = 4;
t=0;
initPoint = [1,0.2];
fun = @(x) bioReactor(t,x,D,X2f);
x = fsolve(fun,initPoint);
disp(x)


%% simulating the non-linear system in part C


% time duration
timeSpan = [0 40];

% equilibrium position
xNought = [0.9943, 1.5141];
% control equilibrium position
D = 0.3;
X2f = 4;

% so here, we're going to simulate this system for varying control inputs
% around the control equilibia, then various perturbations about the
% state equilibria.

% first, the control equilibria
figure
hold on
for i = 0:0.1:1
    [t,y] = ode45(@(t,x) bioReactor(t,x,D+i,X2f+i),timeSpan,xNought);
    plot(t,y)
end
hold off
title('Bio Reactor with Control Perturbations')
xlabel('Time')
ylabel('Output Variables')

% second, the position equilibria
figure
hold on
for i = 0:-0.2:-2
    [t1,y1] = ode45(@(t,x) bioReactor(t,x,D,X2f),timeSpan,xNought+i);
    plot(t1,y1)
end
hold off
title('Bio Reactor with Initial State Perturbations')

%% lastly, open buckle simulation on lineaized state 
% i'm not entirely sure i've been simming the linearized stuff properly in
% the past. So let's have a go at it now.

x0 = [0.2355,-0.9718];
%x0 = [1,1];
D = 0;
A = computeA(xNought,D)
[V,Values] = eig(A)
figure
[t2,y2] = ode45(@(t,z) linearizedBioReactor(t,z,A),timeSpan,x0);
plot(t2,y2)
title('Linearized Bio Reactor')























