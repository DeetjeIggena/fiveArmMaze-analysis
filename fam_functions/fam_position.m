% function returns normalized positions/locations/points 
% depending on maze-version

function [position] = fam_position(position,xmin,xmax,ymin,ymax)

[row,~] = size(position);
% Data-normalization goal position
for r = 1:row
    x(r,1) = datanorm(position(r,1),xmin,xmax);
    y(r,1) = datanorm(position(r,2),ymin,ymax);
end

position = [x, y];
end
