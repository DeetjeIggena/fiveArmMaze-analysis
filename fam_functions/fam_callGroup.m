% function assigns ans returns group-names and group-numbers
% adjust as desired
% depends on ID

function [group,Group] = fam_callGroup(ID, gIDA,gNoA,gNameA,gIDB,gNoB,gNameB)
g = str2double(ID(2));
        % determine group, provide group names as string
    if g == gIDA % group 0=Control, 1=Experimentalgroup
            group = gNoA;
            Group = gNameA;
    elseif g == gIDB
            group = gNoB;
            Group = gNameB;
    else
            group = 100;
            Group = 'F';
            disp('The provided id is out of limits. Group is set to "F"')
    end
end
