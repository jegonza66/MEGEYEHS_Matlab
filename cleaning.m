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