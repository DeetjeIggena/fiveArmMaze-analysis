% function returns the chosen spatial strategy depending on final location

function strategy = fam_strategy(goal, finalZone, condition)

if finalZone == goal && (condition == 3 || condition == 1)
    strategy = 1; % allocentric
elseif (condition == 3 && finalZone == 1) || (condition ==2 && finalZone == 3)
    strategy = 2; % egocentric
else
    strategy = 0;
end
end