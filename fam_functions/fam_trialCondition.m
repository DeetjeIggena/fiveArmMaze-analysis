% function identifies and returns trial-conditions as numbers
% adjust as suits
function [tc] = fam_trialCondition(tIDA, tIDB, tIDC, tIDD, trial_type)

if contains(trial_type,tIDA)
    tc = 0;
elseif contains(trial_type,tIDB)
    tc = 3;
elseif contains(trial_type,tIDC)
    tc = 1;
elseif contains(trial_type,tIDD)
    tc = 2;
else
    tc = 100;
end

end
