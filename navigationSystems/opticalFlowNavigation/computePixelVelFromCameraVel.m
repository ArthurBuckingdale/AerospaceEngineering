function pixelVel = computePixelVelFromCameraVel(x,xData)
% the purpose of this function is to fit a velocity of the camera to match
% that at which the pixels are changing on the image. This will used a 6
% dimensional vector which represents the linear and angular velocity of
% the spacecraft. We have to make some assuptions about the distance to the
% stars. The only important thing to note is that:
% starDist >>>> focal length.

focalLen = 0.893; %meters
focalLength = 0.893 * 1024/0.013312;  % in pixels
tenLightYear = 9.564 * 10.^16; %meters
q=1;
for i = 1:size(xData,1)
    imageJacobianMatrix = [-focalLength/tenLightYear,0,xData(i,1)/tenLightYear,...
        (xData(i,1).*xData(i,2))/focalLength, -(focalLength+((xData(i,1).^2)/focalLength)),xData(i,2);
        0,-focalLength/tenLightYear,xData(i,2)/tenLightYear,(focalLength +((xData(i,2).^2)/focalLength)),...
        -xData(i,1).*xData(i,2)/focalLength,-xData(i,1)];
    pixelVel(q,:) = imageJacobianMatrix * x';
    q=q+1;
end




