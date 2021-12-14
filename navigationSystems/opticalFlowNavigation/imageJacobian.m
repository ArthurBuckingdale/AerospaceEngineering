% the purpose of this is to program the image jacobian matrix in order to
% recover the velocity of the camera from the pixel velocity

focalLength = 0.893; %meters
tenLightYear = 9.461*10^16; %meters
pixelCoords = [500,500];
pixelSpeed = [70/115;-70/115];
cameraSpeed = [0,0,0,0.1,0.1,0.1]';

spaceCraftAngularVeloctity = cross(stateVect(1:3), stateVect(4:6))/norm(stateVect(1:3));

imageJacobianMat = [-focalLength/tenLightYear,0,pixelCoords(1)/tenLightYear,...
    (pixelCoords(1).*pixelCoords(2))/focalLength, -(focalLength+(pixelCoords(1).^2/focalLength)),pixelCoords(2);
    0,-focalLength/tenLightYear,pixelCoords(2)/tenLightYear,(focalLength +(pixelCoords(2).^2/focalLength)),...
    -pixelCoords(1).*pixelCoords(2)/focalLength,-pixelCoords(1)];


pixelSpeed = imageJacobianMat * cameraSpeed;

xData = starCoordinates;
yData = pixelVelocities;
xInitial = [stateVect(4:6),spaceCraftAngularVeloctity];
lb(1:3) = -norm(stateVect(4:6));
lb(4:6) = -norm(spaceCraftAngularVeloctity);
ub(1:3) = norm(stateVect(4:6));
ub(4:6) = norm(spaceCraftAngularVeloctity);

options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');

x = lsqcurvefit(@computePixelVelFromCameraVel,xInitial,xData,yData,lb,ub);