clear all
close all
clc
cd('/mnt/6a6fd40a-e256-4844-8004-0e60d95969e8/MEGEYEHS/Matlab');

% checks user identity and define paths according to local environment
[runpath,code_path,session_path,whoisrunning]=add_paths_matlab_MEG('joaco');

cd(runpath)
fields = fieldnames(code_path);
for ii=1:numel(fields)
    dirtoadd=code_path.(fields{ii}); % this adds a dynamic field name, there's no need to add eeglab 
    addpath(fullfile(dirtoadd),'-end'); % or fieldtrip subdirectories anymore
end

ctfdir                      = session_path.meg;
etdir                       = session_path.et;
bh_dir                      = session_path.behav;
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



fiff_file = 'C:/Users/joaco/OneDrive - The University of Nottingham/MEGEYEHS/Save/Preprocesed_Data/15909001/Subject_15909001_meg.fif';
fiff_file = '/mnt/6a6fd40a-e256-4844-8004-0e60d95969e8/MEGEYEHS/Save/Preprocesed_Data/15909001/Subject_15909001_meg.fif'
cfg = []
cfg.dataset = fiff_file;
data1 = ft_preprocessing(cfg);
ft_datatype(data1)  % returns 'raw'
sessionfilenames = 
evt = mne_read_events(fiff_file);


%%
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
    
    % Load all sessions
    for j=1:length(sessionfilenames)
        % Header
        hdr_name = strcat( 'hdr_',num2str(j));
        eval(sprintf('%s = ft_read_header([ctfdir,filesep,subjname,filesep,sessionfilenames{j}]);',hdr_name))
        
        % Data
        dat_name = strcat( 'dat_',num2str(j));
        eval(sprintf('%s = ft_read_data([ctfdir,filesep,subjname,filesep,sessionfilenames{j}]);',dat_name))
        
        % Events
        evt_name = strcat( 'evt_',num2str(j));
        eval(sprintf('%s = ft_read_event([ctfdir,filesep,subjname,filesep,sessionfilenames{j}]);',evt_name))
        
    end
    
    % Concatenate to 1 variable
    hdr         = hdr_1;
    dat         = dat_1;
    evt         = evt_1;
    nsamples    = hdr_1.nTrials*hdr_1.nSamples; % 1st session events
    
    for j=2:length(sessionfilenames)
        hdr_name = strcat( 'hdr_',num2str(j));
        eval(sprintf('hdr.nTrials = hdr.nTrials + %s.nTrials;', hdr_name))   
        
        dat_name = strcat( 'dat_',num2str(j));
        eval(sprintf('dat = cat(3, dat, %s);', dat_name))
        eval(sprintf('clear %s', dat_name))
        
        evt_name = strcat( 'evt_',num2str(j));

        % shift the sample of the events or triggers in the following block
        hdr_name = strcat( 'hdr_',num2str(j));
        eval(sprintf('for i=1:length(%s); %s(i).sample = %s(i).sample + nsamples; end',evt_name,evt_name,evt_name));
        eval(sprintf('nsamples = nsamples + %s.nTrials*%s.nSamples;',hdr_name,hdr_name)); 
        
        eval(sprintf('evt = cat(1, evt, %s);', evt_name))
        
        eval(sprintf('clear %s',evt_name))
        eval(sprintf('clear %s',hdr_name))
        
    end
    
    clear dat_1 hdr_1 evt_1
    
    % Delete orig data
%     hdr.orig = 'merged';

    
elseif length(sessionfilenames) == 1
    hdr = ft_read_header([ctfdir,filesep,subjname,filesep,sessionfilenames{j}]);
    dat = ft_read_data([ctfdir,filesep,subjname,filesep,sessionfilenames{j}]);
    evt = ft_read_event([ctfdir,filesep,subjname,filesep,sessionfilenames{j}]);
end

% ft_write_data('concatenated.vhdr', dat, 'header', hdr, 'event', evt);

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
xlim([100,1100])
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
blinks_start = find(diff(blinks) == 1);
blinks_end = find(diff(blinks) == -1);
blinks_dur = blinks_end - blinks_start;

% Get actual and fake blinks based on duration condition
actual_blinks = find(blinks_dur > blink_min_samples);
fake_blinks  = find(blinks_dur <= blink_min_samples);

% Define conversion noise intervals around blinks to also fill with nan
start_interval_samples = 5;
end_interval_samples = 12;

plot(ET_MEG.data(1,:), 'DisplayName','Raw')
hold on 

% Fill blinks with nan
for actual_blink_idx=actual_blinks
    blink_interval = blinks_start(actual_blink_idx)-start_interval_samples:blinks_end(actual_blink_idx)+end_interval_samples;
    ET_MEG.data(1, blink_interval) = NaN;
    ET_MEG.data(2, blink_interval) = NaN;
    ET_MEG.data(3, blink_interval) = NaN;
end

plot(ET_MEG.data(1,:), 'g', 'DisplayName','No actual blinks')
hold on

% Remove fake blinks by interpolation
for fake_blink_idx=fake_blinks
    blink_interval = blinks_start(fake_blink_idx)-start_interval_samples:blinks_end(fake_blink_idx)+end_interval_samples;
    interpolation_x = linspace(ET_MEG.data(1, blink_interval(1)), ET_MEG.data(1, blink_interval(end)), length(blink_interval));
    interpolation_y = linspace(ET_MEG.data(2, blink_interval(1)), ET_MEG.data(2, blink_interval(end)), length(blink_interval));
    interpolation_pupil = linspace(ET_MEG.data(3, blink_interval(1)), ET_MEG.data(3, blink_interval(end)), length(blink_interval));
    ET_MEG.data(1, blink_interval) = interpolation_x;
    ET_MEG.data(2, blink_interval) = interpolation_y;
    ET_MEG.data(3, blink_interval) = interpolation_pupil;
    
end

plot(ET_MEG.data(1,:), 'r', 'DisplayName','No blinks')
legend

% Save
if ~exist(fullfile(session_path.preproc_data, subjname))
    mkdir(fullfile(session_path.preproc_data, subjname))
end
save(fullfile(session_path.preproc_data, subjname, 'ET_MEG.mat'),'ET_MEG');

%% 0.3 Define trials based on button box signal

% Button channel name and index
button_label = 'UPPT001';
buttons_idx = find(strcmp(data_MEG.label,button_label));

% Button press start and end
button_diff = diff(data_MEG.trial{1}(buttons_idx,:)) - data_MEG.trial{1}(buttons_idx,1:end-1);

% Get each button press start samples
b1 = find(button_diff == 1);
b2 = find(button_diff == 2);
b4 = find(button_diff == 4);
b8 = find(button_diff == 8);

%% 0.3 bis Define trials based on events

% Defne time array
dat_shape = size(dat);
sfreq = 1/1200;
time = 0:sfreq:dat_shape(2)*dat_shape(3)*sfreq;

%Get events in meg data (only red blue and green)
buttons_idx = find((string({evt.type}) == 'red') | (string({evt.type}) == 'blue') | (string({evt.type}) == 'green') | (string({evt.type}) == 'yellow'));
evt_trial = evt(buttons_idx);
evt_times = time([evt_trial.sample]);
evt_buttons = string({evt_trial.type});

%Check for first trial when first response is not green
for i = 1:length(evt_buttons)
    if evt_buttons(i) ~= 'green'
        first_trial = i;
        break;
    end
end

% Drop events before 1st trial
evt_times = evt_times(first_trial:end);
evt_buttons = evt_buttons(first_trial:end);

% Split events into blocks by green press at begening of block
blocks_start_end = find(evt_buttons == 'green');
blocks_start_end = [0 blocks_start_end length(evt_buttons)]; %-1 is for the first trial 
% Define block starting and ending trial for each block
for i = 1:length(blocks_start_end)-1
    block_bounds(i,:) = [blocks_start_end(i) + 1, blocks_start_end(i+1)];
end

% Load bh data
if length(behavfilenames) > 1
    % get behav file NOT containing 'ORIGINAL'
end

bh_file_path = [bh_dir,filesep,subjname,filesep,behavfilenames{1}];

bh_data = readtable(bh_file_path);

%% borrador
clc
fs=hdr.Fs;

evts    = {evt_1,evt_2,evt_3};
hdrs    = {hdr_1,hdr_2,hdr_3};
s0      = 0;
figure;
    hold on
        for j=1:length(evts)
            evttmp = evts{j};
            hdrtmp = hdrs{j};
        %     EVT_TYPES = unique({evttmp.type});
            EVT_TYPES = {'blue','red'};
            for i = 1:length(EVT_TYPES)
                samples = [evttmp(strcmp({evttmp.type},EVT_TYPES{i})).sample] + s0;
                plot(samples/fs, i*ones(size(samples)),'.')
                fprintf('%s: %d\n',EVT_TYPES{i},length(samples));
            end
            s0 = s0 + hdrtmp.nSamples*hdrtmp.nTrials;
        end
    hold off
    set(gca,'ytick',1:length(EVT_TYPES),'yticklabel',EVT_TYPES)
    set(gca,'ylim',[0 length(EVT_TYPES)+1])
