function [Out] = PolRotate(In, theta, phi)

% Azimuthal rotation angle in radians (in the clockwise direction):
theta = theta*pi/180;

% Elevation rotation angle in radians
phi = phi*pi/180;

% Rotation matrix:
R = [cos(theta) sin(theta)*exp(-1i*phi); -sin(theta)*exp(1i*phi) cos(theta)];

% Applying the rotation to the signals:
Out = (R*In.').';

end