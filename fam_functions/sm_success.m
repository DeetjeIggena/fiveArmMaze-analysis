% function returns whether coordinates are in a certain distance to the
% another position

function [success] = sm_success(distance, threshold)

if distance <= threshold
    success = 1; % yes
else
    success = 0; % no
end
 
end