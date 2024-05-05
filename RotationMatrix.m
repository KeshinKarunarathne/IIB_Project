function [Out] = RotationMatrix(In,TDegree)

% Rotation angle in radians (in the anti-clockwise direction):
Theta = TDegree*pi/180;

% Rotation matrix:
R = [cos(Theta) -sin(Theta); sin(Theta) cos(Theta)];

% Applying the rotation to the signals:
Out(:,1) = R(1,1)*In(:,1) + R(1,2)*In(:,2);
Out(:,2) = R(2,1)*In(:,1) + R(2,2)*In(:,2);

end