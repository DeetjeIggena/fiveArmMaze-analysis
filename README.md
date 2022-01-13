# fiveArmMaze-analysis
Contains script and functions for the analysis of human spatial navigation data, designed for a virtual five armed maze based on the "starmaze"-design (Rondi-reig, et al., 2006, Igloi, et al. 2009)

Used for data analysis in the following manuscripts:


further information:

%% fiveArmMaze-Analysis-script
 @ date 200923 @author Deetje Iggena (deetje.iggena@charite.de)
 @ date 220110 last update
 analysis-script for starmaze version 190819 by Patrizia M Maier
 Matlab R2020b

% Input:
 The analysis-script requires .csv as input-files. Input-files "trackepoint_movment_filename" should contain
 timestamp ("time"), x-position (pos_x), y-position (pos_z), rotation(rot_z)
 
 BE AWARE:
 In case your input is organized differently, please adjust the reading-in process accordingly.

Input-file "trial_results" should contain info about the order of track-files and goal/-task information

% Output:
Output is organized in a structured table and saved as yourTable_currentDate.mat-file
Optional: data is written to a yourTable_currentDate.xlsx

% After starting the script, you are asked to provide the first and the last ID of subjects. The script will search for ID-folders containing the
subfolder 'S001', which should contain the relevant data as single csv.files containing data for single trials. Otherwise you have to change your input depending on your requirements

% After relevant folders have been identified result folder and default-table "fam_currentDate" containing the structured table "fam" will be created.

% In a next step the five armed maze will be created. Herefore, the following data provided by .csv-files is required:
1. minimum as well as maximum x- and y-coordinates e.g. 'fam_minMax.csv'
2. x- and y- coordinaztes for start positions e.g. 'fam_start.csv'
3. x- and y- coordinaztes for goal position e.g. 'fam_goal.csv'
4. x-coordniates for external alleys e.g. 'fam_alleyX.csv'
5. y-coordniates for external alleys e.g. 'fam_alleyY.csv'
6. x-coordniates for inner pentagon e.g. 'fam_pentagonY.csv'
7. y-coordniates for inner pentagon e.g. 'fam_pentagonY.csv'
the five armed maze is saved in matrices and polyshapes

data analysis is embedded in two loops, 1st participant loop, 2nd file loop. Relevant files are identified. Irrelevant files like guidance-trials and log-files are ignored. Trial and group-information is retrieved from trial_results-file and your provided information

% data analysis contains
time-analysis --> latency
trajectory-analysis --> path, distance
final location --> final distance and success
zone-analysis --> coordinates in zones
normalized measures (comparison to ideal trajectory) as error
exploration-analysis --> path-score (number of explored zones)
strategy-analysis --> spatial strategy depending on final chosen location and trial-condition

% Finally, data is saved in yourTable_currentDate.mat
Optional: data is saved as table in yourTable_currentDate.xlsx - file

Quit by entering CTRL+C
