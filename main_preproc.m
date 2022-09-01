clear all
close all
clc

% checks user identity and define paths according to local environment
[runpath,code_path,session_path,whoisrunning]=add_paths_matlab_MEG();

cd(runpath)
fields = fieldnames(code_path);
for ii=1:numel(fields)
    dirtoadd=code_path.(fields{ii}); % this adds a dynamic field name, there's no need to add eeglab 
    addpath(fullfile(dirtoadd),'-end'); % or fieldtrip subdirectories anymore
end

ctfdir                      = session_path.meg;
etdir                       = session_path.et;
matdir                      = session_path.matfiles;

%% Step 0.0 Load data
su = 1;
    subjname            = session_path.subjname{su};
    tmp = dir(fullfile(session_path.meg,subjname,'*.ds')); sessionfilenames = {tmp.name};
    tmp = dir(fullfile(session_path.behav,subjname,'*.csv')); behavfilenames = {tmp.name};

j = 1;
cfg = [];
cfg.matdir      = [matdir,filesep,subjname,filesep,sessionfilenames{j}(1:end-3)];
if ~exist(cfg.matdir,'dir')
    mkdir(cfg.matdir);
end

display('________________________________________________');
display(['subject: ',subjname]);
cfg.myfname         = sessionfilenames{j}([1:(end-3)]);
cfg.dataset         = [ctfdir,filesep,subjname,filesep,sessionfilenames{j}];
cfg.continuous      = 'yes';
cfg.bsfilter        = 'no';
cfg.demean          = 'no';
cfg.detrend         = 'no';
data_MEG            = ft_preprocessing(cfg);

%% Load multiple sessions BETA

if length(sessionfilenames) > 1
    for j=1:length(sessionfilenames)
        var_name = strcat( 'dat_',num2str(j));
        eval(sprintf('%s = ft_read_data([ctfdir,filesep,subjname,filesep,sessionfilenames{j}]);',var_name))
    end
    
    dat = [];
    for j=1:length(sessionfilenames)
        var_name = strcat( 'dat_',num2str(j));
        eval(sprintf('dat = cat(3, dat, %s);', var_name))
        eval(sprintf('clear var_name'))
    end

elseif length(sessionfilenames) == 1
    dat = ft_read_data([ctfdir,filesep,subjname,filesep,sessionfilenames{j}]);
end

%% Step 0.1 Scale ET MEG data

% ET channel names
gaze_x_label = 'UADC001';
gaze_y_label = 'UADC002';
pupil_s_label = 'UADC013';

% Finding the index
numlabel = numel(data_MEG.label); 
for ii=1:numlabel
    if (strcmp(gaze_x_label,data_MEG.label{ii})==1)
        gazex_index = ii;
    end
    if (strcmp(gaze_y_label,data_MEG.label{ii})==1)
        gazey_index = ii;
    end
    if (strcmp(pupil_s_label,data_MEG.label{ii})==1)
        pupil_index = ii;
    end
end

% Plotting the content
figure()
plot(data_MEG.time{1}, data_MEG.trial{1}(gazex_index,:))
hold on;
plot(data_MEG.time{1}, data_MEG.trial{1}(gazey_index,:)-8)
xlim([100,110])
title('MEG ET data')
xlabel('time(s)')
ylabel('gaze')

% Get ET data in MEG
ET_MEG.data = [data_MEG.trial{1}(gazex_index,:); data_MEG.trial{1}(gazey_index,:); data_MEG.trial{1}(pupil_index,:)];
ET_MEG.labels = {'Gaze x', 'Gaze y', 'Pupils'};
ET_MEG.time = data_MEG.time{1};

% Scaling parameters
minvoltage      = -5; %from analog.ini analog_dac_range
maxvoltage      = 5; %from analog.ini analog_dac_range
minrange        = -0.2; %from analog.ini analog_x_range to allow for +/- 20% outside display
maxrange        = 1.2; %from analog.ini analog_x_range to allow for +/- 20% outside display
screenright     = 1919; %OPM lab
screenleft      = 0; 
screentop       = 0;
screenbottom    = 1079; %OPM lab

% Scale
voltageH=ET_MEG.data(find(strcmp(ET_MEG.labels,'Gaze x')),:);
voltageV=ET_MEG.data(find(strcmp(ET_MEG.labels,'Gaze y')),:);
R_h = (voltageH-minvoltage)./(maxvoltage-minvoltage); %voltage range proportion
S_h = R_h.*(maxrange-minrange)+minrange; %proportion of screen width or height
R_v = (voltageV-minvoltage)./(maxvoltage-minvoltage);
S_v = R_v.*(maxrange-minrange)+minrange;
%S_h = ((voltageH-minvoltage)./(maxvoltage-minvoltage)).*(maxrange-minrange)+minrange;
%S_v = ((voltageV-minvoltage)./(maxvoltage-minvoltage)).*(maxrange-minrange)+minrange;
Xgaze = S_h.*(screenright-screenleft+1)+screenleft;
Ygaze = S_v.*(screenbottom-screentop+1)+screentop;    
%now overwrites output in Volts to give the output in pixels
ET_MEG.data = [Xgaze;Ygaze;ET_MEG.data(find(strcmp(ET_MEG.labels,'Pupils')),:)];

% Plotting the scaled content
figure()
plot(ET_MEG.time, ET_MEG.data(find(strcmp(ET_MEG.labels,'Gaze x')),:))
hold on;
plot(ET_MEG.time, ET_MEG.data(find(strcmp(ET_MEG.labels,'Gaze y')),:)-1000)
xlim([100,110])
title('MEG ET scaled data')
xlabel('time(s)')
ylabel('gaze')

% Save
if ~exist(fullfile(session_path.preproc_data, subjname))
    mkdir(fullfile(session_path.preproc_data, subjname))
end
save(fullfile(session_path.preproc_data, subjname, 'ET_MEG.mat'),'ET_MEG');

%% 0.2 Detect Blinks

% Get pupils data to define blinks
meg_pupils_data = ET_MEG.data(find(strcmp(ET_MEG.labels,'Pupils')),:);
blinks = meg_pupils_data < -4.75;

% Minimum blinks duration
blink_min_time = 50; %ms
blink_min_samples = blink_min_time / 1000 * data_MEG.fsample;

% Get blinks start/ end samples and duration
blink_start = find(diff(blinks) == 1);
blink_end = find(diff(blinks) == -1);
blinks_dur = blink_end - blink_start;

% Get actual and fake blinks based on duration condition
actual_blinks = find(blinks_dur > blink_min_samples);
fake_blinks  = find(blinks_dur <= blink_min_samples);

% Define conversion noise intervals around blinks to also fill with nan
interval_samples = 8;

plot(ET_MEG.data(1,:))
hold on 

% Fill blinks with nan
for actual_blinks_idx=actual_blinks
    blink_interval = blink_start(actual_blinks_idx)-interval_samples:blink_end(actual_blinks_idx)+interval_samples;
    ET_MEG.data(1, blink_interval) = NaN;
    ET_MEG.data(2, blink_interval) = NaN;
    ET_MEG.data(3, blink_interval) = NaN;
end

plot(ET_MEG.data(1,:), 'r')

