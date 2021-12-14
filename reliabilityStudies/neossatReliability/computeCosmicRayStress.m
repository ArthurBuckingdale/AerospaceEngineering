function Q = computeCosmicRayStress(shielding, t_sun)
% the purpose of this function is to run the cosmic ray stress code.
for i = 1:7
% now for the latitude. Let's create latitude over the course of one day
time = 1:1:1440;
latitude = 10*i * sind(0.2 * time);

% now for the flux at a given latitude.
cosmicRayFlux = abs(latitude/90)+0.1;

tLat = ones(1,1440);
Q(i,:) = shielding * ((tLat .* cosmicRayFlux)*(t_sun) + 0.1*(tLat .* cosmicRayFlux)*(1 - t_sun));

end
end
