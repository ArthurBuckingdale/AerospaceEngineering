% the purpose of this script is to calculate the rotation matrix given a
% handfull of vectors. It will be estimating the J2000 ECI to the ECEF
% coordinate system.


% j2000
v_a = [-582.70019,2491.37803,6677.801954]'/norm([-582.70019,2491.37803,6677.801954]);
w_a = [-515.6926,2409.0753,6713.46255]'/norm([-582.70019,2491.37803,6677.801954]);
q_a = [-549.2143,2450.3006,6695.8300]'/norm([-549.2143,2450.3006,6695.8300]);

% ECEF
v_b = [1047.21369, 2334.49084, 6677.80195]'/norm([1047.21369, 2334.49084, 6677.80195]);
w_b = [1053.0341, 2227.2632, 6713.4625]'/norm([1053.0341, 2227.2632, 6713.4625]);
q_b = [1050.1836, 2280.948, 6695.8300]'/norm([1050.1836, 2280.948, 6695.8300]);


% first we require a cross product
oneCross = cross(v_a,v_b);

% then we require the dot product
oneDot = dot(v_a,v_b);

% then the cross product matrix of the cross product
crossMat = [0, -oneCross(3), oneCross(2);
            oneCross(3),0,-oneCross(1);
            -oneCross(2),oneCross(1),0];
        
% finally, perform the calculation
rotMat = eye(3) + crossMat + (crossMat^2)*(1/(1+oneDot));

% check the rotation matrix has rotation matrix properties
disp(rotMat')
rotMatInverse = inv(rotMat)
rotMatDet = det(rotMat)

newECEF = rotMat * v_a
dif = v_b - newECEF
























