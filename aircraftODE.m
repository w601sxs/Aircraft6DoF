%% Aircraft ODE
function dx = aircraftODE(t,p) 
dx=zeros(12,1);
    x=p(1); y=p(2); z=p(3);
    vx=p(4); vy=0; vz=p(6);
    psi=p(7); theta=p(8); phi=p(9);
    wx=p(10); wy=p(11); wz=p(12);
%Aircraft parameters   
m=100; g=9.80665;
Ix=1000; Iy=100; Iz=10;
Fx=1000; Fy=0; Fz=0;
Mx=0; My=0; Mz=0;
xcg=10;
xac=12;
xel=17;
yail = 15;
xail=11;
delta_el=5*pi/180;
delta_ail=2*pi/180;


Cel=3; Cw=5; Cd = 1; Cail=2;

% v=[vx - wz*y + wy*z; vy + wz*x - wx*z; vz - wy*x + wx*y];
v=[vx - wz*y + wy*z; 0; vz - wy*x + wx*y];
if(norm(v)==0)
    alpha=0;
else
alpha=asin(-v(3)/norm(v));
end

Lail1= Cail*(alpha - delta_ail);
Lail2= Cail*(alpha + delta_ail);
Lel = -Cel*(alpha - delta_el);
Lw  = Cw * alpha;
Dw = 5 + Cd * alpha;

Mx = Mx + (-yail-0)*Lail1 + (yail-0)*Lail2;
My = My + (xel - xcg)*Lel + (xac - xcg)*Lw + (xail-xcg)*Lail1 + (xail-xcg)*Lail2; 
Fx = Fx + Lw*sin(alpha) + Lel*sin(alpha - delta_el) - Dw*cos(alpha) - Lail1*sin(alpha-delta_ail) - Lail2*sin(alpha+delta_ail) + g*m*sin(theta);
Fy = Fy -g*m*cos(theta)*sin(phi);
Fz = Fz + Lw*cos(alpha) + Lel*cos(alpha - delta_el) + Dw*sin(alpha) + Lail1*cos(alpha-delta_ail) + Lail2*cos(alpha+delta_ail) -g*m*cos(phi)*cos(theta);

Mx= Mx;
My= My;


 
    

    %xd,yd,zd
    dx(1) = vx;
    dx(2) = vy;
    dx(3) = vz;
    
    %xdd, ydd, zdd
    dx(4) = (Fx + m*vy*wz - m*vz*wy)/m - (z*(My - Ix*wx*wz + Iz*wx*wz))/Iy + (Mz*m*y + Ix*m*wx*wy*y - Iy*m*wx*wy*y)/(Iz*m);
    dx(5) = 0;%(Fy - m*vx*wz + m*vz*wx)/m + (z*(Mx + Iy*wy*wz - Iz*wy*wz))/Ix - (Mz*m*x + Ix*m*wx*wy*x - Iy*m*wx*wy*x)/(Iz*m);
    dx(6) = (Fz + m*vx*wy - m*vy*wx)/m - (y*(Mx + Iy*wy*wz - Iz*wy*wz))/Ix + (My*m*x - Ix*m*wx*wz*x + Iz*m*wx*wz*x)/(Iy*m);
    
    %psid, thetad, phid
    dx(7) = (wz*cos(phi) + wy*sin(phi))/cos(theta);
    dx(8) = wy*cos(phi) - wz*sin(phi);
    dx(9) = wx + wz*cos(phi)*tan(theta) + wy*sin(phi)*tan(theta);
    
    %psidd, thetadd, phidd
    dx(10) = (Mx + wy*wz*(Iy - Iz))/Ix;
    dx(11) = (My - wx*wz*(Ix - Iz))/Iy;
    dx(12) = (Mz + wx*wy*(Ix - Iy))/Iz;
    
end

