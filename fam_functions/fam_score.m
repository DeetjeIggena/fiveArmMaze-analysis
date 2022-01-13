% function returns path-score, number of visited zones, max. 20
% last @update 220113

% input: number of entries into relevant areas 
% output: path-score

function score = fam_score(alley_zone_out,alley_zone_in, rectangle_zone, triangle_zone)

[~,col]=size(alley_zone_in);

% determine path score, determine how many zones have been crossed
score = 0;
for c=1:col
    if alley_zone_out(1,c) > 0
        score = score+1;
    end
    if alley_zone_in(1,c) > 0
        score = score+1;
    end
    if rectangle_zone(1,c) > 0
        score = score+1;
    end
    if triangle_zone(1,c) > 0
        score = score+1;
    end
end
end