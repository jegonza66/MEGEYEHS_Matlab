clear all
close all
clc
% cd('/mnt/6a6fd40a-e256-4844-8004-0e60d95969e8/MEGEYEHS/Matlab');
cd('C:\Users\joaco\OneDrive - The University of Nottingham\MEGEYEHS\Matlab')

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
preprocfilenames = [preprocfilenames(end), preprocfilenames(1:end-1)];
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

%% Run visual inspection for every file

% Load fif files
for i=1:length(preprocfilenames)
    cfg  = [];
    cfg.dataset    = fullfile(preproc_dir,filesep,subjname,filesep,preprocfilenames{j});
    cfg.demean = 'yes';    
    data{i,1}  = ft_preprocessing(cfg); 
end

% Concatenate datafiles
cfg = [];
cfg.keepsampleinfo = 'no';
data_planar = ft_appenddata(cfg,data{:});

% Visual inspection to assist rejection
cfg = [];
cfg.viewmode='vertical';
cfg.continuous='yes';
% cfg.ylim=[ -1e-11  1e-11 ];
% cfg.channel = {'MEG01*','MEG05*'};  
ft_databrowser(cfg, data_planar);

% Drop


cfg=[];
cfg.method  = 'summary';
cfg.layout  = 'neuromag306planar.lay';
planar_rjv1  = ft_rejectvisual(cfg, data_planar);

trl_keep=planar_rjv1.trialinfo(1:end,2);
planar_rjv1.trialinfo(:,2)=[];
chan_rej=setdiff(data_planar.label,planar_rjv1.label);
chan_keep=data_planar.label;
chan_keep(find(ismember(data_planar.label(:,1), chan_rej(:,1))))= [];

% Save
save([data_path,'chan_keep.mat'], 'chan_keep')
save([data_path,'trl_keep.mat'], 'trl_keep')
%% ICA

% Downsample
cfg=[]; 
cfg.resamplefs = 200;   
planar_resamp = ft_resampledata(cfg, planar_rjv1);

% Filter
cfg=[];
cfg.hpfreq=1;
cfg.lpfreq=40;
cfg.padding=10.5;
cfg.padtype='zero';
planar_filt=ft_preprocessing(cfg, planar_resamp);

% ICA decomposition
data_rank= rank(planar_filt.trial{1}*planar_filt.trial{1}');
[u,s,v] = svd(planar_filt.trial{1}*planar_filt.trial{1}');
plot(log10(diag(s)),'-r*');

cfg        = [];  
cfg.method = 'runica';  
cfg.numcomponent = data_rank;
planar_comp = ft_componentanalysis(cfg,planar_filt);   

% View components
cfg = [];
cfg.channel = [1:10]; 
cfg.continuous='no';
cfg.viewmode = 'component'; 
cfg.layout = 'neuromag306planar.lay';
ft_databrowser(cfg, planar_comp);

% Reject components
cfg = [];  
cfg.component = [4 37];  
planar_ica = ft_rejectcomponent(cfg, planar_comp, planar_rjv1);

% Plot results
figure
subplot(2,1,1)
plot((1:length(planar_rjv1.trial{1,119}(63,:)))/1000,planar_rjv1.trial{1,119}(63,:),(1:length(planar_ica.trial{1,119}(63,:)))/1000,planar_ica.trial{1,119}(63,:))
xlabel('time (s)')
title('MEG0922 before (blue) and after (red) ICA')
subplot(2,1,2)
plot((1:length(planar_rjv1.trial{1,119}(86,:)))/1000,planar_rjv1.trial{1,119}(86,:),(1:length(planar_ica.trial{1,119}(86,:)))/1000,planar_ica.trial{1,119}(86,:))
xlabel('time (s)')
title('MEG1213 before (blue) and after (red) ICA')

% Re-do visual artifact rejection
planar_ica.trialinfo(1:end,2)=1:length(planar_ica.trialinfo); 

cfg = [];
cfg.viewmode='vertical';
cfg.continuous='no';
artf=ft_databrowser(cfg, planar_ica);

planar_clean=ft_rejectartifact(artf,planar_ica)

trl_keep2=planar_clean.trialinfo(1:end,2);

% Save
save([data_path,'trl_keep2.mat'],'trl_keep2')
save([data_path,'planar_comp.mat'], 'planar_comp')


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
