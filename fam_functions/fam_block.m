% funtion returns the number of the feedback-block depending on block and/
% or trial-number

function block = fam_block(b_num, t_num, fb)
if fb == 0
    if (b_num == 1 && t_num == 10) || (b_num == 2 && t_num == 14) ||...
            (b_num == 3 && t_num == 17) || (b_num == 4 && t_num > 27)
        block = 2;
    else
        block = 1;
    end
else
    block = 0;
end