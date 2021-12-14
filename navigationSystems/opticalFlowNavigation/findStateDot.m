function dot = findStateDot(t,state)
    
    % vecteur state de dimension 6 = [x y z vx vy vz]';
    
    x = state(1);   %(m)
    y = state(2);   %(m)
    z = state(3);   %(m)
    vx = state(4);  %(m/s)
    vy = state(5);  %(m/s)
    vz = state(6);  %(m/s)

    mu = 3.986004418e14; % standart gravitational parameter (m^3/s^2)
    R = 6.37813e6; % Earth radius (m)
    J2 = 1.082626e-3; % dimensional coefficient
    r = norm([x,y,z]); % norm between center of Earth and the satellite position
    
    dot = zeros(4,1); 

    dot(1) = vx;
    dot(2) = vy;
    dot(3) = vz;
    dot(4) = -(mu*x)/r^3 + ((mu*J2*R^2)/2)*(((15*x*z^2)/r^7) - (3*x)/r^5);
    dot(5) = -(mu*y)/r^3 + ((mu*J2*R^2)/2)*(((15*y*z^2)/r^7) - (3*y)/r^5);
    dot(6) = -(mu*z)/r^3 + ((mu*J2*R^2)/2)*(((15*z^3)/r^7) - (9*z)/r^5);
end