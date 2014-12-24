
AircraftDynamics101;

options = odeset('RelTol',1e-10);
[T,Y] = ode45(@aircraftODE,[0 12],[0 0 0 100 0 0 0 0 0 0 0 0],options);


x=Y(:,1); y=Y(:,2); z=Y(:,3);
psi=Y(:,7); theta=Y(:,8); phi=Y(:,9);

P=[y.*(cos(phi).*sin(psi) + cos(psi).*sin(phi).*sin(theta)) + z.*(sin(phi).*sin(psi) - cos(phi).*cos(psi).*sin(theta)) + x.*cos(psi).*cos(theta) ...
  y.*(cos(phi).*cos(psi) - sin(phi).*sin(psi).*sin(theta)) + z.*(cos(psi).*sin(phi) + cos(phi).*sin(psi).*sin(theta)) - x.*cos(theta).*sin(psi) ...
  x.*sin(theta) + z.*cos(phi).*cos(theta) - y.*cos(theta).*sin(phi)];

plot3(P(:,1),P(:,2),P(:,3))
figure, plot(T,Y(:,7),T,Y(:,8),T,Y(:,9)),legend('psi','theta','phi')
