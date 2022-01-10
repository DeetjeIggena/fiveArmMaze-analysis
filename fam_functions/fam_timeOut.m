% function returns whether a time-out occured

function to = fam_timeOut(threshold, time)

if time >= threshold
    to = 1; % yes
else
    to = 0; % no
end