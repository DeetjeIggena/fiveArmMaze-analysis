% function returns the starting-position
function start = fam_trialStart(trial_type)

start1  = 'training'; % alley 1
start11 = 'probe_ego'; % alley 1
start2  = 'probe_allo_startC';
start3  = 'probe_allo_startE';
start4  = 'probe_allo_startG'; %mixed
start5  = 'probe_allo_startI';

if contains(trial_type,start1) || contains(trial_type,start11)
    start = 1;
elseif contains(trial_type,start2)
    start = 2;
elseif contains(trial_type,start3)
    start = 3;
elseif contains(trial_type,start4)
    start = 4;
elseif contains(trial_type,start5)
    start = 5;
end

end