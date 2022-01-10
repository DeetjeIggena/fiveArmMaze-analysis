% function rotates coordniates around the centre to normalize the position
% for further analysis

function [rotatedXY] = fam_rotation(originX,...
    originY, currentX, currentY, x, y, start)
format short
centre_x = 0.5; centre_y = 0.5; % default - normalized data

pointRotation = [x, y];

% calculate angle for rotation matrix 
% determine angle between origin & first start-position 
x1 = [centre_x originX];
y1 = [centre_y originY];
% determine angle between origin & current start-position 
x2 = [centre_x currentX];
y2 = [centre_y currentY];

v1 = [x1(2),y1(2),0] - [x1(1),y1(1),0];
v2 = [x2(2),y2(2),0] - [x1(1),y1(1),0];

angle = round(rad2deg(atan2(norm(cross(v1,v2)), dot(v1,v2))));

if start == 3 || start == 2
    angle = 360-angle;
end

% center of rotation
center = repmat([centre_x, centre_y], 1, 1);
% rotation-matrix
R = [cosd(angle) -sind(angle); sind(angle) cosd(angle)];

point = pointRotation - center;

pointRotated = point*R;

rotatedXY = pointRotated + center;

end
