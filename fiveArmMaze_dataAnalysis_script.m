clear; close all; clc; format compact;

%% fiveArmMaze-Analysis-script
% @ date 200923 @author Deetje Iggena (deetje.iggena@charite.de)
% @ date 220110 last update
% analysis-script for starmaze version 190819
% Matlab R2020b

% Input:
% The analysis-script requires .csv as input-files.
% Input-files "trackepoint_movment_filename" should contain
% timestamp ("time"), x-position (pos_x), y-position (pos_z), rotation
% (rot_z)
% BE AWARE:
% In case your input is organized differently, please adjust
% the reading-in process accordingly.

% Input-file "trial_results" should contain info about the order of
% track-files and goal/-task information

% Output:
% Output is organized in a structured table
% and saved as yourTable_currentDate.mat-file
% Optional: data is written to a yourTable_currentDate.xlsx

% After starting the script, you are asked to provide the first and the
% last ID of subjects. The script will search for ID-folders containing the
% subfolder 'S001', which should contain the relevant data as single csv.files
% containing data for single trials. Otherwise you have to change your
% input depending on your requirements

% After relevant folders have been identified result folder and
% default-table "fam_currentDate" containing the structured table "fam"
% will be created.

% In a next step the five armed maze will be created. Herefore, the
% following data provided by .csv-files is required:
% 1. minimum as well as maximum x- and y-coordinates e.g. 'fam_minMax.csv'
% 2. x- and y- coordinaztes for start positions e.g. 'fam_start.csv'
% 3. x- and y- coordinaztes for goal position e.g. 'fam_goal.csv'
% 4. x-coordniates for external alleys e.g. 'fam_alleyX.csv'
% 5. y-coordniates for external alleys e.g. 'fam_alleyY.csv'
% 6. x-coordniates for inner pentagon e.g. 'fam_pentagonY.csv'
% 7. y-coordniates for inner pentagon e.g. 'fam_pentagonY.csv'
% the five armed maze is saved in matrices and polyshapes

% data analysis is embedded in two loops, 1st group-loops (gN =
% max number of groups), 2nd participant loop (called from group-arrays)
% Relevant files are identified. Irrelevant files like guidance-trials and
% log-files are ignored.
% trial and group-information is retrieved from trial_results-file and your
% provided information

% data analysis contains
% time-analysis --> latency
% trajectory-analysis --> path, distance
% final location --> final distance and success
% zone-analysis --> coordinates in zones
% normalized measures (comparison to ideal trajectory) as error
% exploration-analysis --> path-score (number of explored zones)
% strategy-analysis --> spatial strategy depending on final chosen location
% and trial-condition

% Finally, data is saved in yourTable_currentDate.mat
% Optional: data is saved as table in yourTable_currentDate.xlsx - file

% Quit by entering CTRL+C

%% Provide folder information
currentDirectory      = pwd; % contains data-folder
addpath(genpath(currentDirectory)); % add subfolders containing functions
finalFolderString     = 'S001'; % default --> folder contains starmaze task data

[subjectStart,subjectEnd] = fam_inputSubjectsNo();

% check for existing participants, alternativeley provide participant-list/array
% sort participants into group-arrays

e = 1;

for sub = subjectStart:subjectEnd
    subString = num2str(sub);
    folderIn  = [currentDirectory '\' subString '\' finalFolderString];
    if exist(folderIn)
        partArray(e,1) = sub;
        e = e +1;
    end
end

%% create result & track-folder & table
% construct a path to results-table
formatOut    = 'yymmdd';
date         = datestr(now,formatOut);

resultFolder = [currentDirectory '\fam_results_' date]; %provide the name for your result folder
if ~exist(resultFolder, 'dir')
    mkdir(resultFolder);
    disp('Your folder didn''t exist, a new result-folder was created')
end

% create table for saving data or load existing file
GoalFileName         = ['\fam_' date '.mat'];
GoalFilePath         = fullfile(resultFolder, GoalFileName);
if isfile(GoalFilePath)
    load(GoalFilePath, 'fam');
end

% initialize table of results if not-existing
if ~exist('fam','var')
    fam = [];
    save(GoalFilePath, 'fam');
end

%% Create five-armed maze

% Min-Max-values
values = table2array(readtable('fam_minMax.csv'));
xmin = values(1,1); xmax = values(2,1);
ymin = values(1,2); ymax = values(2,2);

% start-positions(normalized!)
start   = table2array(readtable('fam_start.csv'));
[start] = fam_position(start,xmin,xmax,ymin,ymax);

% goal-positions(normalized!)
goal   = table2array(readtable('fam_goal.csv'));
[goal] = fam_position(goal,xmin,xmax,ymin,ymax);

% coordinates alley-corner (normalized!)
alleyX = table2array(readtable('fam_alleyX.csv'));
[cornerNo,alleyNo] = size(alleyX);
for alley=1:alleyNo
    for corner=1:cornerNo
        alleyX(corner,alley) = datanorm(alleyX(corner,alley),xmin,xmax);
    end
end

alleyY = table2array(readtable('fam_alleyY.csv'));
for alley = 1:alleyNo
    for corner = 1:cornerNo
        alleyY(corner,alley) = datanorm(alleyY(corner,alley),ymin,ymax);
    end
end

% combined pentagon
pentagonX = table2array(readtable('fam_pentagonX.csv'));
pentagonY = table2array(readtable('fam_pentagonY.csv'));

[cP_x,cP_y,cP,pentagonX,pentagonY, inPentagon, outPentagon] = fam_pentagon(alleyX,...
    alleyY,pentagonX, pentagonY,xmin,xmax,ymin,ymax);

% ----------- create polyshapes -------------------------------------------
% alley_polyshape
[alley_full_x,alley_full_y,alley_polyshape, alley_half_out_x, alley_half_out_y, alley_polyshape_1,...
    alley_half_in_x, alley_half_in_y, alley_polyshape_2] = fam_alleyPolyshape(alleyX,alleyY);

% rectangle_polyshape
[rec_x,rec_y,rec] = fam_rectPolyshape(alleyX, alleyY, pentagonX, pentagonY);

% triangle_polyshape
[tri_x,tri_y,tri] = fam_trianglePolyshape(alleyX, alleyY, pentagonX, pentagonY);

% array contains alley & central pentagon polyshapes
polyshape_array = [alley_polyshape_1{1,1} alley_polyshape_2{1,1} alley_polyshape_1{1,2} alley_polyshape_2{1,2}...
    alley_polyshape_1{1,3} alley_polyshape_2{1,3} alley_polyshape_1{1,4} alley_polyshape_2{1,4}...
    alley_polyshape_1{1,5} alley_polyshape_1{1,5} alley_polyshape_2{1,5} cP];

% identify the external goal alley
goalAlley = 0;
for i=1:length(alley_full_x)
    if inpolygon(goal(:,1), goal(:,2), alley_full_x(:,i), alley_full_y(:,i))
        goalAlley = i;
    end
end

n = length(fam) + 1; % counter

participantLength   = length(partArray);

for p = 1: participantLength
    
    % get participant information
    id        = partArray(p);
    idString  = num2str(partArray(p));
    
    % get participant files
    folderIn = [currentDirectory '\' idString '\' finalFolderString];
    
    files = dir(fullfile(folderIn,'*.csv'));
    files = {files.name};
    
    logIndex = find( contains(files,'log'));
    files(:,logIndex)=[];
    
    trialIndex = find( contains(files,'trial_results.csv'));
    files(:,trialIndex)=[];
    
    data_trial=readtable([folderIn, '\trial_results.csv']); % read in trial-file-info
    
    trial = 1;
    
    for f = 1:numel(files)
        
        guide = strcmp(data_trial.trial_type(f,1),'guidance');
        
        if guide == 1
            continue
        end
        
        name = files{f};
        name = fullfile(folderIn, name);
        data = readtable(name);
        
        %% data preparation -> remove pause & first two rows
        % --------------------------------------------------------------------------------------
        if isempty(data)
            trial = trial + 1;
            continue
        end
        
        data(1:2,:) = [];
        data(strcmp(data.gameIsPaused,'True'),:) = [];
        
        if isempty(data)
            trial = trial + 1;
            continue
        end
        
        %% Info depending on single trial: Feedback, Trial-Type/Condition, goal-& startposition
        % --------------------------------------------------------------------------------------
        fam(n).id              = id;
        [fam(n).groupNo,fam(n).groupName]    = fam_callGroup(idString,2,1,'MNE',3,0,'CTR'); % group info
        
        fam(n).trial           = trial;
        fam(n).block_num       = data_trial.block_num(f,1);
        fam(n).trial_num       = data_trial.trial_num(f,1);
        fam(n).trial_in_block  = data_trial.trial_num_in_block(f,1);
        fam(n).trial_type      = data_trial.trial_type(f,1);
        fam(n).trial_condition = fam_trialCondition('training',...
            'startG', 'allo', 'ego', fam(n).trial_type);
        fam(n).feedback        = data_trial.trial_feedback(f,1);
        fam(n).fb              = int8(strcmp(fam(n).feedback,'true'));
        
        fam(n).block_feedback = fam_block(fam(n).block_num,...
            fam(n).trial_num, fam(n).fb);
        
        
        %% select data
        t = data.time;  % time
        x = data.pos_x; % x-coordinate
        y = data.pos_z; % y-coordinate
        
        x = datanorm(x,xmin,xmax); y = datanorm(y,ymin,ymax);   % data normalization for coordinates
        
        dlength = length(x); % get data-length
        
        % Get sampling rate
        % --------------------------------------------------------------------------------------
        samplingRate=zeros((length(data.time)-1),1);
        for i=1:length(data.time)-1
            samplingRate(i)    = data.time(i+1,1) - data.time(i,1);
        end
        avgSamplingRate        = sum(samplingRate)/length(samplingRate);
        
        %% save distinct coordinates
        xStart   = x(1,1);
        yStart   = y(1,1);
        xEnd     = x(end,1);
        yEnd     = y(end,1);
        
        % determine start
        [fam(n).start]   = fam_trialStart(fam(n).trial_type);
        
        [rotatedXY]      = fam_rotation(start(1,1),start(1,2),...
            xStart, yStart, x, y,fam(n).start);
        
        rxEnd   = rotatedXY(end,1); ryEnd = rotatedXY(end,2);
        rxStart = rotatedXY(1,1); ryStart = rotatedXY(1,2);
        
        [rotatedGoal] = fam_rotation(start(1,1),start(1,2),...
            xStart, yStart,goal(1,1),goal(1,2),fam(n).start);
        
        %% data-analysis
        % --------------------------------------------------------------------------------------
        % Time-Analysis using timestamp
        fam(n).time    = t(end,:) - t(1,1); % total amount of time
        
        % get time-out
        fam(n).timeOut = fam_timeOut(88, fam(n).time); % 88 secondsby default
        
        % Coordinate-Analysis using x & y, z
        % -------------------------------------------------------------
        % Path % distance-Analysis
        pathLength          = zeros(1, dlength);
        distGoal            = zeros(1, dlength);
        distFxy             = zeros(1, dlength);
        
        for i=1:dlength-1
            pathLength(i)   = fam_distance(x(i),x(i+1),y(i),y(i+1));% cumulative distance traveled
            distGoal(i)     = fam_distance(x(i),goal(:,1),y(i),goal(:,2)); % cumulative distance to allocentric Goal
            distFxy(i)      = fam_distance(x(i),xEnd,y(i),yEnd); % cumulative distance to allocentric Goal
        end
        
        fam(n).pathLength         = sum(pathLength); % path-length/ distance travelled
        fam(n).avgDistGoal        = sum(distGoal)/dlength; % cumulative distance to allocentric Goal
        fam(n).avgDistFxy         = sum(distFxy)/dlength; % cumulative distance to final xy
        fam(n).velocity           = fam(n).pathLength/fam(n).time;
        
        fam(n).finalDistance      = fam_distance(xEnd,goal(:,1),yEnd,goal(:,2)); % final distance to Goal
        fam(n).success            = fam_success(fam(n).finalDistance, 0.1);
        
        
        %%  Variables depending on starting-positions
        % --------------------------------------------------------------------------------------
        % determine final position --> alley-number or pentagon
        [~,finalAlleyNo,~]    = fam_pointInZone(rxEnd,ryEnd,alley_full_x,alley_full_y,cP_x, cP_y);
        [~,rotGoalAlley,~]    = fam_pointInZone(rotatedGoal(1,1),rotatedGoal(1,2),...
            alley_full_x,alley_full_y,cP_x, cP_y);
        
        % calculate shortest path to goal
        [nodesNo, nodesArray] = fam_shortestPath(start,...
            inPentagon, outPentagon, 1, rotGoalAlley);
        
        % calculate shortest path to final location
        if finalAlleyNo ~= 0
            [~, fxyNodesArray] = fam_shortestPath(start,...
                inPentagon, outPentagon, 1, finalAlleyNo);
            [row,~]= size(fxyNodesArray);
            if row <= 1
                fxyNodesArray = [fxyNodesArray; rxEnd,ryEnd];
            end
            % correction
            fxyNodesArray(end,1) = rxEnd;
            fxyNodesArray(end,2) = ryEnd;
        else
            fxyDist = zeros(length(inPentagon),1);
            for i=1:length(inPentagon)
                fxyDist(i) =   fam_distance(rxEnd,inPentagon(i,1),ryEnd,inPentagon(i,2));
            end
            minDist = min(fxyDist(:));
            [r,~]   = find(fxyDist == minDist);
            r = r + 5;
            [~, fxyNodesArray] = fam_shortestPath(start,...
                inPentagon, outPentagon, 1, r);
            if (r ~= 8 || r ~=9)
                fxyNodesArray(end,1) = rxEnd;
                fxyNodesArray(end,2) = ryEnd;
            elseif r == 8 && (fxyNodesArray(end,1) <= rxEnd)
                fxyNodesArray(end,1) = rxEnd;
                fxyNodesArray(end,2) = ryEnd;
            elseif r == 9 && (fxyNodesArray(end,1) >= rxEnd)
                fxyNodesArray(end,1) = rxEnd;
                fxyNodesArray(end,2) = ryEnd;
            else
                fxyNodesArray = [fxyNodesArray; rxEnd,ryEnd];
            end
        end
        
        % correction
        fxyNodesArray(1,1) = rxStart;fxyNodesArray(1,2) = ryStart;
        nodesArray(1,1) = rxStart; nodesArray(1,2) = ryStart;
        
        
        %% ------- Ideal measures -------------------------------------
        % interpolate data for further analysis
        
        [nodesNum, ~] = size(nodesArray);
        if nodesNum > 1
            [xi_al, yi_al] = fam_interpolation(nodesArray,avgSamplingRate);
            ipLength       = length(xi_al);
            iPathLength    = zeros(1, ipLength);
            iDistGoal      = zeros(1, ipLength);
            
            for i=1:ipLength-1
                iPathLength(i)    = fam_distance(xi_al(i),xi_al(i+1),yi_al(i),yi_al(i+1));% cumulative distance traveled
                iDistGoal(i)      = fam_distance(xi_al(i),rotatedGoal(1,1),yi_al(i),rotatedGoal(1,2)); % cumulative distance to allocentric Goal
            end
            
            fam(n).iPathLength    = sum(iPathLength); % path-length/ distance travelled
            fam(n).iAvgDistGoal   = sum(iDistGoal)/ipLength; % cumulative distance to allocentric Goal
            fam(n).pathError      = fam_error(fam(n).pathLength,fam(n).iPathLength);
            fam(n).distanceError  = fam_error(fam(n).avgDistGoal,fam(n).iAvgDistGoal);
        else
            if fam(n).pathLength         <= 0.1
                fam(n).iPathLength        = 0;
                fam(n).iAvgDistGoal       = 0;
                fam(n).pathError          = 0;
                fam(n).distanceError      = 0;
            else
                fam(n).iPathLength     = 0.1; % movement-tolerance
                fam(n).iAvgDistGoal    = 0.05; %  movement-tolerance
                fam(n).pathError       = fam_error(fam(n).pathLength,fam(n).iPathLength);
                fam(n).distanceError   = fam_error(fam(n).avgDistGoal,fam(n).iAvgDistGoal);
            end
        end
        
        [fxyNodesNum, ~] = size(fxyNodesArray);
        z = fxyNodesArray(1,:) ~= fxyNodesArray(2,:);
        if fxyNodesNum > 2 && z(1,1)==1 && z(1,2)==1
            [fxi_al, fyi_al]    = fam_interpolation(fxyNodesArray,avgSamplingRate);
            ipLength            = length(fxi_al);
            fiPathLength        = zeros(1, ipLength);
            fiDistLoc           = zeros(1, ipLength);
            
            for i=1:ipLength-1
                fiPathLength(i)    = fam_distance(fxi_al(i),fxi_al(i+1),fyi_al(i),fyi_al(i+1));% cumulative distance traveled
                fiDistLoc(i)       = fam_distance(fxi_al(i),rxEnd,fyi_al(i),ryEnd); % cumulative distance to allocentric Goal
            end
            
            fam(n).fiPathLength    = sum(fiPathLength); % path-length/ distance travelled
            fam(n).fiAvgDistGoal   = sum(fiDistLoc)/ipLength; % cumulative distance to allocentric Goal
            fam(n).fxyPathError    = fam_error(fam(n).pathLength,fam(n).fiPathLength);
            fam(n).fxyDistanceError   = fam_error(fam(n).avgDistGoal,fam(n).fiAvgDistGoal);
        else
            if fam(n).pathLength          <= 0.1
                fam(n).fiPathLength        = 0;
                fam(n).fiAvgDistGoal       = 0;
                fam(n).fxyPathError        = 0;
                fam(n).fxyDistanceError    = 0;
            else
                fam(n).fiPathLength        = 0.1; % movement-tolerance
                fam(n).fiAvgDistGoal       = 0.05; %  movement-tolerance
                fam(n).fxyPathError        = fam_error(fam(n).pathLength,fam(n).fiPathLength);
                fam(n).fxyDistanceError    = fam_error(fam(n).avgDistGoal,fam(n).fiAvgDistGoal);
            end
        end
        %% Zone-Analysis
        [zone_alley_abs,zone_alley_rel,zone_alley_time,...
            zone_alley_entry] = fam_coordinatesArea(x,y,...
            alley_full_x,alley_full_y,fam(n).time);
        
        zone_alley_abs_all    = sum(zone_alley_abs);
        zone_alley_rel_all    = sum(zone_alley_rel);
        zone_alley_time_all   = sum(zone_alley_time);
        zone_alley_entry_all  = sum(zone_alley_entry);
        
        [zone_rectangle_abs, ~,~, ~]  = fam_coordinatesArea(x,y,...
            rec_x,rec_y, fam(n).time);
        [zone_triangle_abs, ~,~, ~]   = fam_coordinatesArea(x,y,...
            tri_x,tri_y, fam(n).time);
        [zone_pentagon_abs,~,~, ~]    = fam_coordinatesArea(x,y,...
            tri_x,tri_y, fam(n).time);
        
        [zone_alleyOut_abs, ~, ~, ~]  = fam_coordinatesArea(x,y,...
            alley_half_out_x,alley_half_out_y, fam(n).time);
        [zone_alleyIn_abs, ~, ~, ~]   = fam_coordinatesArea(x,y,...
            alley_half_in_x,alley_half_in_y, fam(n).time);
        
        [goalAlleyRel, goalAlleyTime, goalAlleyEntry] = fam_zoneCorrection(goalAlley, fam(n).start, ...
            zone_alley_rel, zone_alley_time, zone_alley_entry);
        
        fam(n).zoneAlleyEntryCorrect   = zone_alley_entry(fam(n).start) + goalAlleyEntry;
        fam(n).zoneAlleyTimeCorrect    = zone_alley_time(fam(n).start) + goalAlleyTime;
        fam(n).zoneAlleyRelCorrect     = zone_alley_rel(fam(n).start) + goalAlleyRel;
        
        fam(n).zoneAlleyRelIncorrect   = zone_alley_rel_all - fam(n).zoneAlleyRelCorrect;
        fam(n).zoneAlleyTimeIncorrect  = zone_alley_time_all - fam(n).zoneAlleyTimeCorrect;
        fam(n).zoneAlleyEntryIncorrect = zone_alley_entry_all - fam(n).zoneAlleyEntryCorrect;
        
        %% Exploration
        [~,fam(n).finalAlleyNo,~] = fam_pointInZone(xEnd,yEnd,...
            alley_full_x,alley_full_y,cP_x, cP_y);
        
        fam(n).spatialStrategy    = fam_strategy(goalAlley,...
            fam(n).finalAlleyNo, fam(n).trial_condition);
        
        % path-score/ exploration of zones
        fam(n).pathScore          = fam_score(zone_alleyOut_abs,...
            zone_alleyIn_abs,...
            zone_rectangle_abs,zone_triangle_abs);
        
        trial = trial + 1;
        n = n+1;
    end
end
% end

% % -----------------------------------------------------------------------
% % Save data
% % -----------------------------------------------------------------------

save(GoalFilePath, 'fam', '-append');

%%
% % -----------------------------------------------------------------------
% % Save data xlsx
% % -----------------------------------------------------------------------

fam_writeToXLSX(fam, resultFolder, 'fam')

close all
clear variables
