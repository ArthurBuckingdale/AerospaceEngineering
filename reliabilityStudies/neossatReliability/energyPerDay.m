function EperDay = energyPerDay(area,solarAngle,distanceFromSun,secondsInSunlightPerDay,albedoes)
%the purpose of this function is to compute the energy per day gained by a
%spacecraft. 
solarLuminosity = 3.828*10^26;
EperDay = area*cos(solarAngle)*(1-albedoes)*...
    (solarLuminosity/(4*pi*(distanceFromSun.^2)))*secondsInSunlightPerDay;
end
