% the purpose of this script is to create a base model for the temperature
% of a spacecraft. 
clear all
close all

albedoes = 0.1:0.1:0.9;
area = 0.1:0.1:0.9;
solarAngle = 0:0.1:pi/2;
secondsInSunlightPerDay = 1800:1800:86400;
distanceFromSun = 30*10^9:5*10^9:300*10^9;

solarLuminosity = 3.828*10^26;


tic
for i = 1:length(area)
    for j = 1:length(solarAngle)
        for h = 1:length(albedoes)
            for k = 1:length(distanceFromSun)
                for l = 1:length(secondsInSunlightPerDay)
                    energyInSpacecraft(i,j,h,k,l) = area(i)*cos(solarAngle(j))*(1-albedoes(h))*...
                        (solarLuminosity/(4*pi*(distanceFromSun(k).^2)))*secondsInSunlightPerDay(l);
                end
            end
        end
    end
end
time = toc;

figure
imagesc(energyInSpacecraft(:,:,end,end,end))
title('Heat Map of Phase Angle(rad*10) and Area(m^{2}*10)')
xlabel('Phase Angle (rad*10)')
ylabel('Area (m^{2}*10)')

figure
imagesc(reshape(energyInSpacecraft(:,end,:,end,end),length(area),[]))
title('Heat Map of albedo(*10) and Area(*10)')
xlabel('Albedo (*10)')
ylabel('Area (m^{2}*10)')

figure
plot(distanceFromSun,reshape(energyInSpacecraft(end,end,end,:,end),length(distanceFromSun),1))
title('Total Energy as a function of distance from sun')
ylabel('Energy(J)')
xlabel('Distance from Sun (m)')

figure
plot(reshape(energyInSpacecraft(end-3,end-3,end-3,:,:),length(distanceFromSun),[]))
title('Heat Map of Seconds in Sun per day And Distance From Sun')
xlabel('Distance from Sun(30*10^9 m)')
ylabel('Energy (J)')
legend('Varying Values of Seconds in Sun Per Day')


specificHeat = 500;
SCMass = 72; %kg
%EperDay = energyPerDay(area,solarAngle,distanceFromSun,secondsInSunlightPerDay,albedoes)
EperDay = energyPerDay(1.2,0,150*10^9,84600,0.7);
radiationCooling = 350; %w
radiatedPerDay = 84600*radiationCooling;
temp = (EperDay-radiatedPerDay)/(specificHeat*SCMass); %K





         