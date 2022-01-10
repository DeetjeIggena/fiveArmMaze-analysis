% functions returns error/deviation/accuracy, actual value deviates from 
% ideal value in percent
% minus values are corrected
function error = fam_error(actual,ideal)

if actual == 0
    error = 0;
else
    error =((actual-ideal)/ideal)*100;
    if error <= 0
        error = error*(-1);
    end
end

end