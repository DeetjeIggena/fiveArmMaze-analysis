clear; close all; clc; format compact;

%% Motor control task analysis script
% @ date 200923 @author Deetje Iggena (deetje.iggena@charite.de)
% @ date 220113 last update
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
% default-table "mct_currentDate" containing the structured table "mct"
% will be created.

% data analysis is embedded in two loops, 1st participant-loop,
% 2nd file loop (.csv files containing the relevant data)
% Relevant files are identified. Irrelevant files like guidance trials,
% trial results and log-files are ignored.
% trial and group-information is retrieved from your
% provided information

% data analysis contains
% time-analysis --> latency
% trajectory-analysis --> path, distance, velocity

% Finally, data is saved in yourTable_currentDate.mat
% Optional: data is saved as table in yourTable_currentDate.xlsx - file

% Quit by entering CTRL+C

%% Provide folder information
currentDirectory      = pwd; % contains data-folder
addpath(genpath(currentDirectory)); % add subfolders containing functions
finalFolderString     = 'S002'; % default --> folder contains starmaze task data

%% select participants
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

resultFolder = [currentDirectory '\mct_results' date]; %provide the name for your result folder
if ~exist(resultFolder, 'dir')
    mkdir(resultFolder);
    disp('Your folder didn''t exist, a new result-folder was created')
end


% create table for saving data
% construct the Goal-file/ load existing Goal-table
GoalFileName         = '\mct.mat';
GoalFilePath         = fullfile(resultFolder, GoalFileName);
if isfile(GoalFilePath)
    load(GoalFilePath, 'mct');
end
% initialize table of results if not-existing
if ~exist('mct','var')
    mct = [];
    save(GoalFilePath, 'mct');
end

n = length(mct) +1;

%n = 1;
participantLength   = length(partArray);

for p = 1: participantLength
    
    % get participant information
    id                          = partArray(p);
    idString                    = num2str(partArray(p));
    [groupNo,groupName]         = fam_callGroup(idString,2,1,'MNE',3,0,'CTR'); % group info
    
    % get participant files
    folderIn    = [currentDirectory '\' idString '\' finalFolderString];
    
    files       = dir(fullfile(folderIn,'*.csv'));
    files       = {files.name};
    
    logIndex    = find( contains(files,'log'));
    files(:,logIndex)   = [];
    
    trialIndex          = find( contains(files,'trial_results.csv'));
    files(:,trialIndex) = [];
    
    data_trial          = readtable([folderIn, '\trial_results.csv']); % read in trial-file-info
    
    trial = 1;
    
    for f = 1:numel(files)
        
        % data preparation
        name = files{f};
        name = fullfile(folderIn, name);
        data = readtable(name);
        
        % remove pause & first two rows
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
        
        % get participant and trial info
        mct(n).id          = id;
        mct(n).groupNo     = groupNo;
        mct(n).groupName   = groupName;
        mct(n).trialNum    = trial;
        
        % get data
        t = data.time;  % time
        x = data.pos_x; % x-coordinate
        y = data.pos_z; % y-coordinate
        
        dlength = length(x); % get data-length
        
        %% data-analysis
        % --------------------------------------------------------------------------------------
        % Time-Analysis using timestamp
        mct(n).time = t(end,:) - t(1,1); % total amount of time
        
        % Coordinate-Analysis using x & y, z
        pathLength       = zeros(1, dlength);
        
        for i=1:dlength-1
            pathLength(i)   = fam_distance(x(i),x(i+1),y(i),y(i+1));% cumulative distance traveled
        end
        
        mct(n).pathLength   = sum(pathLength); % path-length/ distance travelled
        mct(n).velocity     = mct(n).pathLength/mct(n).time;
        
        trial = trial + 1;
        n     = n +1;
    end
end


%%  Save data
% % ---------------------------------------------------------------------------------------------------------------------------
save(GoalFilePath, 'mct', '-append');


%%  Save data xlsx
% % ---------------------------------------------------------------------------------------------------------------------------
fam_writeToXLSX(mct, resultFolder, 'mct')

close all
clear variables