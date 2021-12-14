% the purpose of this function is to create a model which estimates the
% stress on a space based platform for cosmic rays. In the LATEX document
% you can find out how this model is built. 
close all
clear all

% initialising some constant parameters.
shielding = 0.5;
t_sun = 0.9;

n=1;
for i = [0.2,0.4,0.6,0.9]
    Q(:,:,n) = computeCosmicRayStress(shielding, i);
    n=n+1;
end

time = 1:1:length(Q(1,:,1));

figure
plot(time, Q(:,:,1))
title('Stress Factor as a Function of Time')
xlabel('Time (min)')
ylabel('Stress')
legend('10 deg','20 deg','30 deg','40 deg','50 deg','60 deg','70 deg')

figure
plot(time, reshape(Q(1,:,:),[],4))
title('Stress as a function of time in sun as percentage of orbit')
xlabel('Time (min)')
ylabel('Stress')
legend('0.2','0.4','0.6','0.9')



















