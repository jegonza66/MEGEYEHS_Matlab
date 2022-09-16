%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% preprocess_hybridsearch_ET_behav.m
%% Script to analyse hybrid_search experiment
%% MJI & JG, version 20.05.22
%% MJI, version April 13, 2022
%% MJI, version 26.02.2020
%% Dependencies: add_paths_matlab_hybridsearch, behavioural files, ET files
%% Output: This code generates the following MATLAB structures:
% trial: Structure with behavioural information on each trial (210) and
% some summary performance on trial.behav
% subject: Structure with participant ID and basic demographics
% psypy: Structure with Psychopy specific information as well as all the
% triggers from psychopy in substructure psypy.timing
% Exp: Structure with participant ID and session locations

clear
close all
clc

extractEyeData      = 1; %perform EM analysis
saveFiles           = 1; %1=saves analysed data
overwrite_behav     = 1; %1=forces overwrite of behavioural data generation
[runpath,codes_path,session_path,whoisrunning]=add_paths_matlab_MEG;
%checks user identity and define paths according to local environment
if ~exist(runpath,'dir'); mkdir(runpath); end
cd(runpath)

cd(runpath)
fields = fieldnames(codes_path);
for ii=1:numel(fields)
    dirtoadd=codes_path.(fields{ii});        % this adds a dynamic field name
                                            % there's no need to add eeglab
                                            % or fieldtrip subdirectories anymore
    addpath(fullfile(dirtoadd),'-end');
end

Nsubj = length(session_path.subjname);
clear trial,clear subject;clear psypy;
%% Step 0. Behavior 
%% Step 0.1. Check if I have all the files, and complete the ones that are left
% JG making variables for all subjects plots
subjects_accu_all = [];
subjects_accu_tpres = [];

subjects_accu_all_n1 = [];
subjects_accu_all_n2 = [];
subjects_accu_all_n4 = [];

subjects_accu_tpres_n1 = [];
subjects_accu_tpres_n2 = [];
subjects_accu_tpres_n4 = [];

subjects_RT_Tpres_n1 = [];
subjects_RT_Tpres_n2 = [];
subjects_RT_Tpres_n4 = [];

subjects_RT_all_n1 = [];
subjects_RT_all_n2 = [];
subjects_RT_all_n4 = [];

for su = 1:Nsubj %subj=4 (S104) no ET data
    %%
    Exp = [];
    Exp.subjname            = session_path.subjname{su};
    tmp = dir(fullfile(session_path.behav, Exp.subjname,'*.csv')); behavfilenames = {tmp.name};
    %check for exact name of behavioural file
    dir_content=dir([fullfile(session_path.behav, Exp.subjname,behavfilenames{su})]);
    Exp.sessionfilenames    = fullfile(dir_content.folder,dir_content.name);    
    display('________________________________________________');
    display(['subject: ',Exp.subjname]);
    Exp.procdir      = fullfile(session_path.out,Exp.subjname);
    if ~exist(Exp.procdir,'dir'); mkdir(Exp.procdir); end
      
    %check if behavioral files already exist
    outfile = fullfile(Exp.procdir,['trial_data_S',Exp.subjname]);
    if (exist(fullfile([outfile '.mat']),'file') && overwrite_behav==0)
        fprintf('%s: Behavioral data file found\n',Exp.subjname)
        localExp=Exp;
        load(fullfile(outfile));
        %now overwrite fields that could be different
        Exp = localExp;
    else
        fprintf('%s: ERROR: Behavioral data cannot be found. Generating it\n',Exp.subjname)
        % Behaviour
        if su < 3
            % Run old behavioural analysis for first to participants who
            % used old version of psychopy expreiment
            [trial,subject,psypy,accu_all_trials_mean,accu_Tpres_trials_mean,accu_all_n1,accu_all_n2,accu_all_n4,accu_Tpres_n1,accu_Tpres_n2,accu_Tpres_n4,RT_Tpres_n1,RT_Tpres_n2,RT_Tpres_n4,RT_all_trials_n1,RT_all_trials_n2,RT_all_trials_n4] = ana_behav_hybridsearch(Exp,1);
        else
            [trial,subject,psypy,accu_all_trials_mean,accu_Tpres_trials_mean,accu_all_n1,accu_all_n2,accu_all_n4,accu_Tpres_n1,accu_Tpres_n2,accu_Tpres_n4,RT_Tpres_n1,RT_Tpres_n2,RT_Tpres_n4,RT_all_trials_n1,RT_all_trials_n2,RT_all_trials_n4] = ana_behav_hybridsearch_NEW(Exp,1);
        end
        subjects_accu_all(su) = accu_all_trials_mean;
        subjects_accu_tpres(su) = accu_Tpres_trials_mean;
        
        subjects_accu_all_n1(su) = accu_all_n1;
        subjects_accu_all_n2(su) = accu_all_n2;
        subjects_accu_all_n4(su) = accu_all_n4;
        
        subjects_accu_Tpres_n1(su) = accu_Tpres_n1;
        subjects_accu_Tpres_n2(su) = accu_Tpres_n2;
        subjects_accu_Tpres_n4(su) = accu_Tpres_n4;
        
        subjects_RT_Tpres_n1(su) = RT_Tpres_n1;
        subjects_RT_Tpres_n2(su) = RT_Tpres_n2;
        subjects_RT_Tpres_n4(su) = RT_Tpres_n4;
        
        subjects_RT_all_n1(su) = RT_all_trials_n1;
        subjects_RT_all_n2(su) = RT_all_trials_n2;
        subjects_RT_all_n4(su) = RT_all_trials_n4;
        
        %Read stim positions and identities from externally generated
        %annotations
        [trial]       = merge_stim_positions(session_path,trial);     
        if saveFiles
            save(outfile,'trial','subject','psypy','Exp', 'session_path');
        end
    end

    %% process eye tracking data
    etdir=fullfile(session_path.et);
    etdir_out=fullfile(session_path.out,session_path.subjname{su});
    if extractEyeData
        et_file = fullfile(etdir,Exp.subjname,['et_data_',Exp.subjname]);
        
        % Raw ET
        if exist(fullfile([et_file '.edf']),'file')
            fprintf('%s: The raw EDF file exists\n',Exp.subjname)
        else
            fprintf('%s: ERROR: The raw EDF cannot be found\n',Exp.subjname)
        end
        %Ascii
        if ~exist(fullfile([et_file '.asc']))
            fprintf('%s: Generating ASC file\n',Exp.subjname)
            cd(etdir)
            eval(['!edf2asc '  Exp.subjname '/et_data_' Exp.subjname '.edf'])
            cd(runpath)
        end
        %% Load Eye Data      
        if ~exist(fullfile(etdir_out,['et_data_S',Exp.subjname,'_all_eyelink.mat']))            
            [all] = fun_extract_all_eye_data([et_file,'.asc'], etdir); 
            %[mall ma]  =  read_eyetracker_data_ASCII_edit()
            save(fullfile(etdir_out,['et_data_S',Exp.subjname,'_all_eyelink.mat']), 'all')
        else
            load(fullfile(etdir_out,['et_data_S',Exp.subjname,'_all_eyelink.mat']))
        end
   
        ET.synctime             = all.msgtime(cellfun(@(x) any(strfind(x,'ETSYNC')),all.msg));
        ET.eyemtime_start       = all.msgtime(cellfun(@(x) any(strfind(x,'ETSYNC 152')),all.msg));
        %ET.eyemtime_start_msg   = 'EEG trigger ETSYNC 150';

        ET.eyemtime_stop        = all.msgtime(cellfun(@(x) any(strfind(x,'ETSYNC 151')),all.msg));
        
        ET.vstime_start         = all.msgtime(cellfun(@(x) any(strfind(x,'ETSYNC 250')),all.msg));
        ET.vstime_start_msg     = 'EEG trigger ETSYNC 255';
        ET.vstime_stop          = all.msgtime(cellfun(@(x) any(strfind(x,'ETSYNC 251')),all.msg));
        
        ET.memtime_start        = all.msgtime(cellfun(@(x) any(strfind(x,'ETSYNC 200')),all.msg));
        ET.memtime_stop         = all.msgtime(cellfun(@(x) any(strfind(x,'ETSYNC 201')),all.msg));
        
        ET.fix_start            = all.msgtime(cellfun(@(x) any(strfind(x,'ETSYNC 50')),all.msg));
        ET.p_port_ini_exp       = all.msgtime(cellfun(@(x) any(strfind(x,'ETSYNC 100')),all.msg));
        ET.p_port_ini_exp_msg   = 'EEG trigger ETSYNC 100';

        %emfile = fullfile(procdir,Exp.subjname,[sessionfilenames '_all_eyelink.mat']);
        %% Parse VS trials
        t_bgn = ET.vstime_start;
        t_end = ET.vstime_stop;
        ET.eye          = all.ojo;
        ET.srate        = all.srate;
        %dbstop if error
        ET.VS.eyedata = fun_parse_eyedata(all,t_bgn,t_end);    
        %trial_to_plot='first';
        trial_to_plot = 'positions_annotated';
        Exp.screenXpixels = 1920; %S10, S11 %BE CAREFUL S12
        Exp.screenYpixels = 1080; %S10, S11 %BE CAREFUL S12
        %cfg.screenXpixels = 1980; %S12
        %cfg.screenYpixels = 1080; %S12
        Ntrials = trial.Ntrials;
        switch trial_to_plot
            case 'random'
                trials_to_plot =  randi(trial.Ntrials);
            case 'all'
                trials_to_plot = 1:trial.Ntrials;
            case 'first'
                trials_to_plot = 1;
            case 'positions_annotated'
                trials_to_plot=[];
                for itt=1:Ntrials
                    if(~isempty(trial.item{itt}))
                        trials_to_plot = [trials_to_plot itt];
                    end
                end
            otherwise
                trials_to_plot = 13;
                warning('trial_to_plot not defined. Plotting trial 10')
        end
        %S105
        %trials_to_plot=[13 22 35 48 53 55 76 83 96 112 114 121 125 130 142 166 168 172 178 197 203 205]; %:210; %1:20
        plot_scanpath(ET,trial,Exp,session_path,trials_to_plot)
        if saveFiles
            outfile = fullfile(Exp.procdir,['ET_data_S',Exp.subjname]);
            save(outfile,'ET', 'Exp', 'session_path')
        end
    end
end


% end
% 
% % Plot all subjects Acc vs Rt and fit
% subjects_accu_all_trials = cat(1,subjects_accu_all_n1, subjects_accu_all_n2, subjects_accu_all_n4);
% subjects_RT_all_trials = cat(1,subjects_RT_all_n1, subjects_RT_all_n2, subjects_RT_all_n4);
% 
% P = polyfit(subjects_RT_all_trials,subjects_accu_all_trials,1);
% % yfit = P(1)*x+P(2);
% 
% yfit = polyval(P, subjects_RT_all_trials);          % Estimated  Regression Line
% SStot = sum((subjects_accu_all_trials-mean(subjects_accu_all_trials)).^2);                    % Total Sum-Of-Squares
% SSres = sum((subjects_accu_all_trials-yfit).^2);                       % Residual Sum-Of-Squares
% Rsq = 1-SSres/SStot;   
% 
% 
% f = figure(); clf
% f.Position = [100 100 1000 400]; 
% plot(subjects_RT_all_trials, subjects_accu_all_trials, 'bo', subjects_RT_all_trials, yfit,'r-')
% ylim([0,1.5])
% xlabel('RT')
% ylabel('Accuracy')
% mylabel = sprintf('Linear Fit: r2=%.2f',Rsq)
% h_legend=legend('Data', mylabel);
% set(h_legend,'FontSize',13,'Location','SouthWest');
% mytitle = sprintf('Subjects Accuracy vs RT All');
% title(mytitle,'FontSize',17);
% outfile=fullfile(session_path.out_ana,'Scatter_Subjects_Accuracy_vs_RT');
% print(outfile,'-dpng','-r0')
% 
% 
% % Acc per MSS ALL
% box_plot_accu_data_all = transpose([subjects_accu_all_n1;subjects_accu_all_n2;subjects_accu_all_n4]);
% f = figure(); clf
% boxplot(box_plot_accu_data_all,'Labels',{'N=1','N=2','N=4'})
% ylabel('Accuracy')
% mytitle = sprintf('Subjects Accuracy per MSS All Trials');
% title(mytitle,'FontSize',17);
% outfile=fullfile(session_path.out_ana,'BoxPlot_Subjects_Accuracy_MSS_All');
% print(outfile,'-dpng','-r0')
% 
% 
% % RT per MSS ALL
% box_plot_RT_data_all = transpose([subjects_RT_all_n1;subjects_RT_all_n2;subjects_RT_all_n4]);
% f = figure(); clf
% boxplot(box_plot_RT_data_all,'Labels',{'N=1','N=2','N=4'})
% ylabel('RT')
% mytitle = sprintf('Subjects RT per MSS All Trials');
% title(mytitle,'FontSize',17);
% outfile=fullfile(session_path.out_ana,'BoxPlot_Subjects_RT_MSS_All');
% print(outfile,'-dpng','-r0')
% 
% 
% % Acc per MSS Tpres
% box_plot_accu_data_Tpres = transpose([subjects_accu_Tpres_n1;subjects_accu_Tpres_n2;subjects_accu_Tpres_n4]);
% f = figure(); clf
% boxplot(box_plot_accu_data_Tpres,'Labels',{'N=1','N=2','N=4'})
% ylabel('Accuracy')
% mytitle = sprintf('Subjects Accuracy per MSS Tpres Trials');
% title(mytitle,'FontSize',17);
% outfile=fullfile(session_path.out_ana,'BoxPlot_Subjects_Accuracy_MSS_Tpres');
% print(outfile,'-dpng','-r0')
% 
% 
% % RT per MSS Tpres
% box_plot_RT_data_Tpres = transpose([subjects_RT_Tpres_n1;subjects_RT_Tpres_n2;subjects_RT_Tpres_n4]);
% f = figure(); clf
% boxplot(box_plot_RT_data_Tpres,'Labels',{'N=1','N=2','N=4'})
% ylabel('RT')
% mytitle = sprintf('Subjects RT per MSS Tpres Trials');
% title(mytitle,'FontSize',17);
% outfile=fullfile(session_path.out_ana,'BoxPlot_Subjects_RT_MSS_Tpres');
% print(outfile,'-dpng','-r0')
% 
% 
% Plot all subjects Acc vs Rt
% f = figure(); clf
% f.Position = [100 100 1000 400];
% plot(subjects_RT_all_n4, subjects_accu_all_n4, 'o')
% ylim([0,1])
% xlabel('RT')
% ylabel('Accuracy')
% mytitle = sprintf('Subjects Accuracy vs RT for All trials of N=4 MSS');
% title(mytitle,'FontSize',17);
% outfile=fullfile(session_path.out_ana,'Scatter_Subjects_Accuracy_vs_RT_All_N4');
% print(outfile,'-dpng','-r0')
% 
% 
% % Plot all subjects mean accuraccy
% f = figure(); clf
% f.Position = [100 100 1000 400]; 
% plot(1:Nsubj, subjects_accu_all, 'o',1:Nsubj, subjects_accu_tpres,'*')
% % % plot(1:Nsubj, subjects_accu_tpres, 'o')
% ylim([0,1])
% h_legend=legend('All Trials', 'Tpres Trials');
% set(h_legend,'FontSize',13,'Location','SouthWest');
% xticks([1:Nsubj])
% xticklabels([session_path.subjname])
% mytitle = sprintf('Subjects Accuracy');
% title(mytitle,'FontSize',17);
% outfile=fullfile(session_path.out_ana,'Scatter_Subjects Accuracy');
% print(outfile,'-dpng','-r0')  
% 
% 
% % Plot all subjects accuraccy for N=4 Memory set size
% f = figure(); clf
% f.Position = [100 100 1000 400]; 
% plot(1:Nsubj, subjects_accu_all_n4, 'o', 1:Nsubj, subjects_accu_tpres_n4,'*')
% % % plot(1:Nsubj, subjects_accu_tpres, 'o')
% ylim([0,1])
% h_legend=legend('All Trials', 'Tpres Trials');
% set(h_legend,'FontSize',13,'Location','SouthWest');
% xticks([1:Nsubj])
% xticklabels([session_path.subjname])
% mytitle = sprintf('Subjects Accuracy N=4 MSS');
% title(mytitle,'FontSize',17);
% outfile=fullfile(session_path.out_ana,'Scatter_Subjects Accuracy_N4_MSS');
% print(outfile,'-dpng','-r0')  
% 
% 
% % Plot all subjects RT for N=2 and N=4 MSS
% f = figure(); clf
% f.Position = [100 100 1000 400]; 
% plot(1:Nsubj, subjects_RT_Tpres_n2, 'o', 1:Nsubj, subjects_RT_Tpres_n4,'*')
% % % plot(1:Nsubj, subjects_accu_tpres, 'o')
% ylim([0,10])
% h_legend=legend('N=2 MSS', 'N=4 MSS');
% set(h_legend,'FontSize',13,'Location','SouthWest');
% xticks([1:Nsubj])
% xticklabels([session_path.subjname])
% mytitle = sprintf('Subjects RT N=2 vs N=4 MSS');
% title(mytitle,'FontSize',17);
% outfile=fullfile(session_path.out_ana,'Scatter_Subjects_RT_N2_vs_N4_MSS');
% print(outfile,'-dpng','-r0')  
