% I will copy and organice myself the previous MEG analysis codes.
% (2022-05-09, JG) I started from a copy of preprocess_MEG_wPreexp_7_multiuser.m

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load, filter, accomodate and merge multiple files in one dataset per participant
%% File path
% from super_preprocess_simple.m
clear all
close all
clc

% checks user identity and define paths according to local environment
[runpath,code_path,session_path,whoisrunning]=add_paths_matlab_MEG();

cd(runpath)
fields = fieldnames(code_path);
for ii=1:numel(fields)
    dirtoadd=code_path.(fields{ii});        % this adds a dynamic field name
                                            % there's no need to add eeglab
                                            % or fieldtrip subdirectories anymore
    addpath(fullfile(dirtoadd),'-end');
end

ctfdir                      = session_path.raw;
etdir                       = session_path.rawet;
matdir                      = session_path.matfiles;

%Empty subjects are removed and info copied to session_path.exclude  
% [session_path] = remove_empty_subjects(session_path);

%% Step 0.1 Define data and Check the triggers channel (UPPT002) 
ID = [];

su = 1 %1:length(session_path.subjname)
subjname            = session_path.subjname{su};
subjcode            = session_path.subjcode{su};
sessionfilenames    = session_path.sessionfilenames{su};
behavfilenames      = session_path.behavfilenames{su};

j = 1
ID = [ID;[su j]];        

cfg = [];
cfg.matdir      = [matdir,filesep,subjname,filesep,sessionfilenames{j}(1:end-3)];
if ~exist(cfg.matdir,'dir')
    mkdir(cfg.matdir);
end

display('________________________________________________');
display(['subject: ',subjname]);
cfg.myfname         = sessionfilenames{j}([1:(end-3)]);
cfg.dataset         = [ctfdir,filesep,sessionfilenames{j}];
cfg.continuous      = 'yes';
%             cfg.channel         = {'UPPT002'}; % (JK, 21/03/2018) I include the two eye tracking channels (UADC013 and UADC014)
                                   % (JK, 28/08/2018) The reference is important for proper work of 3rd gradient correction
cfg.bsfilter        = 'no';
cfg.demean          = 'no';
cfg.detrend         = 'no';
data_MEG                = ft_preprocessing(cfg);

%% Checking content of individual Eye channels
gaze_x_label = 'UADC001';
gaze_y_label = 'UADC002';
pupil_s_label = 'UADC013';
button_1_label = 'UPPT001';
button_2_label = 'UPPT002';

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

%% Load edf data
etdir=fullfile(session_path.rawet);
etdir_out=fullfile(session_path.out,session_path.subjname{su});
et_file = fullfile(etdir,Exp.subjname,['et_data_',subjname]);

% Raw ET
if exist(fullfile([et_file '.edf']),'file')
    fprintf('%s: The raw EDF file exists\n',subjname)
else
    fprintf('%s: ERROR: The raw EDF cannot be found\n',subjname)
end

%Ascii
if ~exist(fullfile([et_file '.asc']))
    fprintf('%s: Generating ASC file\n',subjname)
    cd(etdir)
    eval(['!edf2asc '  subjname '/et_data_' subjname '.edf'])
    cd(runpath)
end

if ~exist(fullfile(etdir_out,['et_data_S',subjname,'_all_eyelink.mat']))            
    [all] = fun_extract_all_eye_data([et_file,'.asc'], etdir); 
    save(fullfile(etdir_out,['et_data_S',subjname,'_all_eyelink.mat']), 'all')
else
    load(fullfile(etdir_out,['et_data_S',subjname,'_all_eyelink.mat']))
end

[new_all] = read_eyetracker_data_ASCII_edit([et_file,'.asc'], etdir);

et_time = all.samples(:,1);
et_gaze_x = all.samples(:,2);
et_gaze_y = all.samples(:,3);
et_pupil_s = all.samples(:,4); 


figure()
plot(et_time, et_gaze_x)
title('ET data')
xlabel('time(s)')
ylabel('gaze x')

%% joint plot

time_dif = 1308700;

norm_MEG_gaze_x = data_MEG.trial{1}(gazex_index,:) - min(data_MEG.trial{1}(gazex_index,:));
norm_MEG_gaze_x = norm_MEG_gaze_x/max(norm_MEG_gaze_x);

norm_et_gaze_x = et_gaze_x - min(et_gaze_x);
norm_et_gaze_x = norm_et_gaze_x / max(norm_et_gaze_x);

figure()
plot(data_MEG.time{1}, norm_MEG_gaze_x)
hold on;
plot(data_MEG.time{1} + 105.2, norm_et_gaze_x(1:length(data_MEG.time{1}))*3 -1)
h_legend=legend('MEG data', 'ET data');
set(h_legend,'FontSize',10,'Location','SouthWest');
ylim([-1,1])
xlim([100,200])
