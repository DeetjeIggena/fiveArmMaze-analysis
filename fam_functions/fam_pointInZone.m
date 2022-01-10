% function returns about final position in alley or pentagon
% input information about the relevant areas and their x-y-coordinates

function [final_alley, final_alley_no,...
    final_pentagon] = fam_pointInZone(x,y,area_x,area_y,cP_x, cP_y)

[~,col]     = size(area_x);
final_alley = zeros(1,col);

final_alley_no = 0;
for c = 1:col
    final_alley(c) = inpolygon(x,y,area_x(:,c), area_y(:,c));
    if inpolygon(x,y,area_x(:,c), area_y(:,c))
        final_alley_no = c ;
    end
    final_pentagon = inpolygon(x,y,cP_x, cP_y);
    
end