function phi = findJacobian(state)

    % vecteur state de dimension 6 = [x y z vx vy vz]';
    
    x = state(1);   %(m)
    y = state(2);   %(m)
    z = state(3);   %(m)
    
    mu = 3.986004418e14; % standart gravitational parameter (m^3/s^2)
    R = 6.37813e6; % Earth radius (m)
    J2 = 1.082626e-3; % dimensional coefficient

    %ax
    daxdx = 3*mu*J2*R^2*(4*x^4+3*x^2*y^2-27*x^2*z^2-y^4+3*y^2*z^2+4*z^4)/(2*(x^2+y^2+z^2)^(9/2));
    daxdy = 15*mu*x*y*R^2*J2*(x^2+y^2-6*z^2)/(2*(x^2+y^2+z^2)^(9/2));
    daxdz = 45*z*(x^2+y^2-(4/3)*z^2)*x*mu*J2*R^2/(2*(x^2+y^2+z^2)^(9/2));

    %ay
    daydx = 15*mu*x*y*R^2*J2*(x^2+y^2-6*z^2)/(2*(x^2+y^2+z^2)^(9/2));
    daydy = -3*mu*J2*R^2*(x^4-3*x^2*y^2-3*x^2*z^2-4*y^4+27*y^2*z^2-4*z^4)/(2*(x^2+y^2+z^2)^(9/2));
    daydz = 45*z*(x^2+y^2-(4/3)*z^2)*y*mu*J2*R^2/(2*(x^2+y^2+z^2)^(9/2));

    %az
    dazdx = 45*z*(x^2+y^2-(4/3)*z^2)*x*mu*J2*R^2/(2*(x^2+y^2+z^2)^(9/2));
    dazdy = 45*z*(x^2+y^2-(4/3)*z^2)*y*mu*J2*R^2/(2*(x^2+y^2+z^2)^(9/2));
    dazdz = -9*(x^4+(2*y^2-8*z^2)*x^2 + y^4-8*y^2*z^2+(8/3)*z^4)*mu*J2*R^2/(2*(x^2+y^2+z^2)^(9/2));


    
    phi = [0 0 0 1 0 0;
           0 0 0 0 1 0;
           0 0 0 0 0 1;
           daxdx daxdy daxdz 0 0 0;
           daydx daydy daydz 0 0 0;
           dazdx dazdy dazdz 0 0 0;];

end