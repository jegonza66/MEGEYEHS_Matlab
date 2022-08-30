clear all
close all
clc

% checks user identity and define paths according to local environment
[runpath,code_path,session_path,whoisrunning]=add_paths_matlab_MEG('joaco');

cd(runpath)
fields = fieldnames(code_path);
for ii=1:numel(fields)
    dirtoadd=code_path.(fields{ii});        % this adds a dynamic field name
                                            % there's no need to add eeglab
                                            % or fieldtrip subdirectories anymore
    addpath(fullfile(dirtoadd),'-end');
end

ctfdir                      = session_path.meg;
etdir                       = session_path.et;
matdir                      = session_path.matfiles;

%% Step 0.0 Load data
su = 1;
    subjname            = session_path.subjname{su};
    subjcode            = session_path.subjcode{su};
    tmp = dir(fullfile(session_path.meg,session_path.subjname{su},'*.ds')); sessionfilenames = {tmp.name}; %session_path.sessionfilenames{su};
    tmp = dir(fullfile(session_path.behav,session_path.subjname{su},'*.csv')); behavfilenames = {tmp.name}; %session_path.behavfilenames{su};

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
%             cfg.channel         = {'UPPT002'}; % (JK, 21/03/2018) I include the two eye tracking channels (UADC013 and UADC014)
                                   % (JK, 28/08/2018) The reference is important for proper work of 3rd gradient correction
cfg.bsfilter        = 'no';
cfg.demean          = 'no';
cfg.detrend         = 'no';
data_MEG                = ft_preprocessing(cfg);

%% Step 0.1 Scale ET MEG data

% ET channel names
gaze_x_label = 'UADC001';
gaze_y_label = 'UADC002';
pupil_s_label = 'UADC013';


% finding the index
numlabel = numel(data_MEG.label); %422
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


%plotting the content
figure()
plot(data_MEG.time{1}, data_MEG.trial{1}(gazex_index,:))
hold on;
plot(data_MEG.time{1}, data_MEG.trial{1}(gazey_index,:)-8)
xlim([100,110])
title('MEG ET data')
xlabel('time(s)')
ylabel('gaze x')

% Get ET data in MEG
ET_MEG.data = [data_MEG.trial{1}(gazex_index,:); data_MEG.trial{1}(gazey_index,:); data_MEG.trial{1}(pupil_index,:)];
ET_MEG.labels = {'Gaze x', 'Gaze y', 'Pupils'};
ET_MEG.time = data_MEG.time{1};

% Scale
minvoltage      = -5; %from analog.ini analog_dac_range
maxvoltage      = 5; %from analog.ini analog_dac_range
minrange        = -0.2; %from analog.ini analog_x_range to allow for +/- 20% outside display
maxrange        = 1.2; %from analog.ini analog_x_range to allow for +/- 20% outside display
screenright     = 1919; %OPM lab
screenleft      = 0; 
screentop       = 0;
screenbottom    = 1079; %OPM lab

voltageH=ET_MEG.data(find(strcmp(ET_MEG.labels,'Gaze y')),:);
voltageV=ET_MEG.data(find(strcmp(ET_MEG.labels,'Gaze x')),:);
eye_diam_meg=ET_MEG.data(find(strcmp(ET_MEG.labels,'Pupils')),:);
R_h = (voltageH-minvoltage)./(maxvoltage-minvoltage); %voltage range proportion
S_h = R_h.*(maxrange-minrange)+minrange; %proportion of screen width or height
R_v = (voltageV-minvoltage)./(maxvoltage-minvoltage);
S_v = R_v.*(maxrange-minrange)+minrange;
%S_h = ((voltageH-minvoltage)./(maxvoltage-minvoltage)).*(maxrange-minrange)+minrange;
%S_v = ((voltageV-minvoltage)./(maxvoltage-minvoltage)).*(maxrange-minrange)+minrange;
Xgaze = S_h.*(screenright-screenleft+1)+screenleft;
Ygaze = S_v.*(screenbottom-screentop+1)+screentop;    
%now overwrites output in Volts to give the output in pixels
ET_MEG.data = [Xgaze;Ygaze;eye_diam_meg];


% Plotting the scaled content
figure()
plot(ET_MEG.time, ET_MEG.data(find(strcmp(ET_MEG.labels,'Gaze x')),:))
hold on;
plot(ET_MEG.time, ET_MEG.data(find(strcmp(ET_MEG.labels,'Gaze y')),:)-8)
xlim([100,110])
title('MEG ET scaled data')
xlabel('time(s)')
ylabel('gaze x')

%% 
