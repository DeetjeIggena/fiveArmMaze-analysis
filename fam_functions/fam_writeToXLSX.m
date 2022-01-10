% function writes table data to xlsx.file

function fam_writeToXLSX(fam, resultFolder)

formatOut = 'yymmdd';
date      = datestr(now,formatOut);

% name of excel-file
file_name = ['fam_' date '.xlsx'];
new_file = fullfile(resultFolder, file_name);

S = [fam(:)];
writetable(struct2table(S),new_file);


