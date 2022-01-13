% function returns corrected values for rel. amount of coordinates, time
% and entries into zones

function [rel, time, entry] = fam_zoneCorrection(goal, start, zoneRel, zoneTime, zoneEntry)

if goal == start
    entry      = 0;
    time       = 0;
    rel        = 0;
else
    entry      = zoneEntry(goalAlley);
    time       = zoneTime(goalAlley);
    rel        = zoneRel(goalAlley);
end

end