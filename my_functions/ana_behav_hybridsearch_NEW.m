function [trial,subject,psypy,accu_all_trials_mean,accu_Tpres_trials_mean,accu_all_n1,accu_all_n2,accu_all_n4,accu_Tpres_n1,accu_Tpres_n2,accu_Tpres_n4,RT_Tpres_n1,RT_Tpres_n2,RT_Tpres_n4,RT_all_trials_n1,RT_all_trials_n2,RT_all_trials_n4] = ana_behav_hybridsearch(Exp,doplots)
%ANA_BEHAV_HYBRIDSEARCH analyses behaviour of hybrid search experiment
%MJI, version 19.02.2020

if(nargin == 1)
    doplots = 0;
else
    doplots = 1;
end

psypy   = [];
subject = [];
trial = [];
%% Analyse behaviour (works much better with xlsx file: CHANGE IN PSYCHOPY)

[num,txt,raw] = xlsread(Exp.sessionfilenames);

%[num,txt,raw] = xlsread(Exp.sessionfilenames,'trials');
%[status,sheets, xlFormat] = xlsfinfo(Exp.sessionfilenames)
total_length  = 210+6*7+3; %this should be Ntrials + 6 x Nblocks + 3
headers       = txt(1,:);
trial_idx_txt = [8:37, 44:73, 80:109, 116:145, 152:181, 188:217, 224:253];
trial_idx_num = trial_idx_txt - 1;

% vs_filenames
lookup = 'searchimage';
rows = find(cellfun(@(c) ischar(c) && strcmp(c, lookup), headers));
trial.vs_filenames = {txt{trial_idx_txt,rows}};
% T_filenames
lookup = 'st5';
rows = find(cellfun(@(c) ischar(c) && strcmp(c, lookup), headers));
trial.T_filenames  = {txt{trial_idx_txt,rows}};
% D1
lookup = 'st1';
rows = find(cellfun(@(c) ischar(c) && strcmp(c, lookup), headers));
trial.D1_filenames  = {txt{trial_idx_txt,rows}};
% D2
lookup = 'st2';
rows = find(cellfun(@(c) ischar(c) && strcmp(c, lookup), headers));
trial.D2_filenames  = {txt{trial_idx_txt,rows}};
% D3
lookup = 'st3';
rows = find(cellfun(@(c) ischar(c) && strcmp(c, lookup), headers));
trial.D3_filenames  = {txt{trial_idx_txt,rows}};
% D4
lookup = 'st4';
rows = find(cellfun(@(c) ischar(c) && strcmp(c, lookup), headers));
trial.D4_filenames  = {txt{trial_idx_txt,rows}};
% Stim present (cell 5 x Ntotal)
lookup = 'st1_cat';
rows = find(cellfun(@(c) ischar(c) && strcmp(c, lookup), headers));
stim_present{1}        =  cell2mat(txt(trial_idx_txt,rows))=='D'; % Changed to 3 because empty row 2. Check why before it worked before with filenames
lookup = 'st2_cat';
rows = find(cellfun(@(c) ischar(c) && strcmp(c, lookup), headers));
stim_present{2}        =  cell2mat(txt(trial_idx_txt,rows))=='D';
lookup = 'st3_cat';
rows = find(cellfun(@(c) ischar(c) && strcmp(c, lookup), headers));
stim_present{3}        =  cell2mat(txt(trial_idx_txt,rows))=='D';
lookup = 'st4_cat';
rows = find(cellfun(@(c) ischar(c) && strcmp(c, lookup), headers));
stim_present{4}        =  cell2mat(txt(trial_idx_txt,rows))=='D';
lookup = 'st5_cat';
rows = find(cellfun(@(c) ischar(c) && strcmp(c, lookup), headers));
stim_present{5}        =  cell2mat(txt(trial_idx_txt,rows))=='T';

trial.stim_present = stim_present;

trial.Ntrials          = numel(trial.vs_filenames);

lookup = 'Nstim';
rows = find(cellfun(@(c) ischar(c) && strcmp(c, lookup), headers));
trial.Nstim        = num(trial_idx_num,rows-2);

% stimulus duration
lookup = 'stDur';
rows = find(cellfun(@(c) ischar(c) && strcmp(c, lookup), headers));
trial.stDurMem     = num(trial_idx_num,rows-2); %this stim duration depends on the number of stimuli
trial.stDurVS      = 10; %VS time search can last up to 10 secs

% responses
lookup = 'key_resp.corr';
rows = find(cellfun(@(c) ischar(c) && strcmp(c, lookup), headers));
trial.respcorr     = num(trial_idx_num,rows-2);

lookup = 'key_resp.rt';
rows = find(cellfun(@(c) ischar(c) && strcmp(c, lookup), headers));
trial.rt           = num(trial_idx_num,rows-2);

trial.behav.perfo = nanmean(trial.respcorr); %mean perfo
trial.behav.rt_mean = nanmean(trial.rt); %mean rt

trial.ind_setsize1=find(trial.Nstim==1);
trial.ind_setsize2=find(trial.Nstim==2);
trial.ind_setsize4=find(trial.Nstim==4);

trial.set_size=[1 2 4];
trial.behav.accu_setsize = [nanmean(trial.respcorr(trial.ind_setsize1)) nanmean(trial.respcorr(trial.ind_setsize2))...
    nanmean(trial.respcorr(trial.ind_setsize4))];
trial.behav.rt_setsize = [nanmean(trial.rt(trial.ind_setsize1)) nanmean(trial.rt(trial.ind_setsize2)) ...
    nanmean(trial.rt(trial.ind_setsize4))];

Tpres=stim_present{5};

trial.behav.ind_setsize1_Tpres=find(trial.Nstim==1 & Tpres==1);
trial.behav.ind_setsize2_Tpres=find(trial.Nstim==2 & Tpres==1);
trial.behav.ind_setsize4_Tpres=find(trial.Nstim==4 & Tpres==1);

trial.behav.accu_setsize_Tpres = [nanmean(trial.respcorr(trial.behav.ind_setsize1_Tpres)) nanmean(trial.respcorr(trial.behav.ind_setsize2_Tpres)) nanmean(trial.respcorr(trial.behav.ind_setsize4_Tpres))];
trial.behav.rt_setsize_Tpres = [nanmean(trial.rt(trial.behav.ind_setsize1_Tpres)) nanmean(trial.rt(trial.behav.ind_setsize2_Tpres)) nanmean(trial.rt(trial.behav.ind_setsize4_Tpres))];

%% Reading table directly to access information from incomplete columns
T=readtable(Exp.sessionfilenames); %this is needed to get the info from incomplete columns

psypy.timing.tracker_message_started        = T.tracker_message_started(2);
%psypy.timing.p_port_ini_exp_started         = T.p_port_ini_exp_started(1);
%psypy.timing.p_port_ini_exp_started_value   = 100; %this is the EEG trigger value at the beg of exp
%psypy.timing.p_port_end_experiment_started  = T.p_port_end_experiment_started(end);
%psypy.timing.p_port_end_experiment_started_value = 255; %this is the EEG trigger value at end of exp
%psypy.timing.p_port_started                 = T.p_port_started(6:6:end);
psypy.timing.T_search_img_2_started           = T.search_img_2_started(7:7:end); 
%psypy.timing.p_port_ini_eyemap_started      = T.p_port_ini_eyemap_started; %this is a tricky variable.
%First five numbers are eyemap stimuli, followed by a NaN. Then it signals passing through
%eyemap look without stopping (close to beginning of memory trial),
%psypy.timing.p_port_ini_eyemap_started_value = 150; 
psypy.psychopyVersion                       = T.psychopyVersion{1};
psypy.frameRate                             = T.frameRate(1);

subject.participant= T.participant(1);
subject.gender= T.gender{1};
subject.age= T.age(1);

%%
%%
%figure;

if doplots
    figure(11); clf
    clc
    subplot(2,1,1)
    hData=scatter(trial.set_size,trial.behav.accu_setsize);
    accu_all_trials_mean = mean(trial.behav.accu_setsize);
    accu_all_n1 = trial.behav.accu_setsize(1);
    accu_all_n2 = trial.behav.accu_setsize(2);
    accu_all_n4 = trial.behav.accu_setsize(3);
    lsline;
    hold on
    hData2=scatter(trial.set_size,trial.behav.accu_setsize_Tpres);
    accu_Tpres_trials_mean = mean(trial.behav.accu_setsize_Tpres);
    accu_Tpres_n1 = trial.behav.accu_setsize_Tpres(1);
    accu_Tpres_n2 = trial.behav.accu_setsize_Tpres(2);
    accu_Tpres_n4 = trial.behav.accu_setsize_Tpres(3);
    lsline;
    set(hData                         , ...
        'Marker'          , 'o'         , ...
        'SizeData'      , 60           , ...
        'MarkerEdgeColor' , 'none'      , ...
        'MarkerFaceColor' , [.75 .75 1] );
    set(hData2                         , ...
        'Marker'          , 's'         , ...
        'SizeData'      , 60           , ...
        'MarkerEdgeColor' , 'none'      , ...
        'MarkerFaceColor' , [1 0 0] );
    xlabel('Memory Set Size','FontSize',14,'FontWeight','bold')
    ylabel('Accuracy','FontSize',14,'FontWeight','bold')
    set(gcf,'color','w')
    All_trials_legend = sprintf('All trials (mean:%.2f)', accu_all_trials_mean);
    Tpres_trials_legend = sprintf('T pres (mean:%.2f)', accu_Tpres_trials_mean);
    h_legend=legend(All_trials_legend, 'linear fit',Tpres_trials_legend);
    set(h_legend,'FontSize',10);
    legend('boxoff')
    set(h_legend,'Location','SouthWest');
    ylim([0 1])
    %legend('FontSize',18)
    box off;
    %
    subplot(2,1,2)
    hData=scatter(trial.set_size,trial.behav.rt_setsize);
    RT_all_trials_n1 = trial.behav.rt_setsize(1);
    RT_all_trials_n2 = trial.behav.rt_setsize(2);
    RT_all_trials_n4 = trial.behav.rt_setsize(3);
    RT_Tpres_n1 = trial.behav.rt_setsize(1);
    RT_Tpres_n2 = trial.behav.rt_setsize(2);
    RT_Tpres_n4 = trial.behav.rt_setsize(3);
    lsline;
    hold on
    hData2=scatter(trial.set_size,trial.behav.rt_setsize_Tpres);
    lsline;
    set(hData                         , ...
        'Marker'          , 'o'         , ...
        'SizeData'      , 60           , ...
        'MarkerEdgeColor' , 'none'      , ...
        'MarkerFaceColor' , [.75 .75 1] );
    set(hData2                         , ...
        'Marker'          , 's'         , ...
        'SizeData'      , 60           , ...
        'MarkerEdgeColor' , 'none'      , ...
        'MarkerFaceColor' , [1 0 0] );
    b=polyfit(trial.set_size,exp(trial.behav.rt_setsize_Tpres),1);  % ln y = mx + b
    yhat=log(polyval(b,trial.set_size));       % compute fitted result
    hLyhat=plot(trial.set_size,yhat,'r--');        % add to the plot; s
    %set(gca,'yscale','log')
    xlabel('Memory Set Size','FontSize',14,'FontWeight','bold')
    ylabel('RT (s)','FontSize',14,'FontWeight','bold')
    set(gcf,'color','w')
    box off;
    ylim([0 7])
    currentFigure = gcf;
    mytitle = sprintf('Behavioral data, subjID:%s',Exp.subjname);
    title(currentFigure.Children(end), mytitle,'FontSize',20);
    %
    %%
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    outfile=fullfile(Exp.procdir,['scatter_behav_subj_',Exp.subjname]);
    print(outfile,'-dpng','-r0')
    
end
end

