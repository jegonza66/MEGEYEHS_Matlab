clear all
close all
clc
% cd('/mnt/6a6fd40a-e256-4844-8004-0e60d95969e8/MEGEYEHS/Matlab');
% cd('C:\Users\joaco\OneDrive - The University of Nottingham\MEGEYEHS\Matlab')

% checks user identity and define paths according to local environment
[runpath,code_path,session_path,whoisrunning]=add_paths_matlab_MEG();

cd(runpath)
fields = fieldnames(code_path);
for ii=1:numel(fields)
    dirtoadd=code_path.(fields{ii}); % this adds a dynamic field name, there's no need to add eeglab 
    addpath(fullfile(dirtoadd),'-end'); % or fieldtrip subdirectories anymore
end

preproc_dir                 = session_path.preproc;
ctfdir                      = session_path.meg;
etdir                       = session_path.et;
bh_dir                      = session_path.behav;
matdir                      = session_path.matfiles;

%% 0 Load data
%% Step 0.1 Define paths and data files
su = 1;
    subjname            = session_path.subjname{su};
    tmp = dir(fullfile(session_path.preproc,subjname,'*meg*.fif')); preprocfilenames = {tmp.name};
    preprocfilenames = [preprocfilenames(end), preprocfilenames(1:end-1)]
    tmp = dir(fullfile(session_path.preproc,subjname,'*eve*.fif')); eventfilenames = {tmp.name};
    tmp = dir(fullfile(session_path.preproc,subjname,'*eve_map.csv')); event_map_filename = {tmp.name};
    tmp = dir(fullfile(session_path.meg,subjname,'*.ds')); sessionfilenames = {tmp.name};
    tmp = dir(fullfile(session_path.behav,subjname,'*.csv')); behavfilenames = {tmp.name};

j = 1;
cfg = [];
cfg.matdir      = [matdir,filesep,subjname,filesep,sessionfilenames{j}(1:end-3)];
if ~exist(cfg.matdir,'dir')
    mkdir(cfg.matdir);
end

%% Step 0.2 Load multiple sessions BETA

if length(preprocfilenames) > 1
    
    % Load all sessions
    for j=1:length(preprocfilenames)
        % Header
        hdr_name = strcat( 'hdr_',num2str(j));
        eval(sprintf('%s = ft_read_header([preproc_dir,filesep,subjname,filesep,preprocfilenames{j}]);',hdr_name))
        
        % Data
        dat_name = strcat( 'dat_',num2str(j));
        eval(sprintf('%s = ft_read_data([preproc_dir,filesep,subjname,filesep,preprocfilenames{j}]);',dat_name))   
    end
    
    % Concatenate to 1 variable
    hdr         = hdr_1;
    dat         = dat_1;
    
    for j=2:length(preprocfilenames)
        hdr_name = strcat( 'hdr_',num2str(j));
        eval(sprintf('hdr.nTrials = hdr.nTrials + %s.nTrials;', hdr_name))   
        
        dat_name = strcat( 'dat_',num2str(j));
        eval(sprintf('dat = cat(2, dat, %s);', dat_name))
        eval(sprintf('clear %s', dat_name)) 
    end
    
clear dat_1 hdr_1
    
elseif length(preprocfilenames) == 1
    hdr = ft_read_header([preproc_dir,filesep,subjname,filesep,preprocfilenames{j}]);
    dat = ft_read_data([ctfdir,filesep,subjname,filesep,preprocfilenames{j}]);
end

if length(eventfilenames) > 1

    % Load all sessions
    for j=1:length(eventfilenames)      
        % Events
        evt_name = strcat( 'evt_',num2str(j));
        eval(sprintf('%s = ft_read_event([preproc_dir,filesep,subjname,filesep,preprocfilenames{j}]);',evt_name))
    end
    
    % Concatenate to 1 variable
    evt         = evt_1;
    
    for j=2:length(eventfilenames)
        evt_name = strcat( 'evt_',num2str(j));
    end
    
    clear dat_1 hdr_1 evt_1
    
elseif length(eventfilenames) == 1
    evt = mne_read_events([preproc_dir,filesep,subjname,filesep,eventfilenames{1}]);
end

% Map of events ids and definition
event_map = readmatrix([preproc_dir,filesep,subjname,filesep,event_map_filename{1}]);