function [trial]       = merge_stim_positions(session_path,trial)
%merge_stim_positions Read stimuli positions from Matlab file, find
%matching trials and incorporate positions to trial structure
%MJI, version 19.02.2020

load(session_path.positions_file); %output is a structure called 'pos'
%loop over all pos elements
%this needs to be revised...template matching?
for ii=1:trial.Ntrials
    %pos_red = pos( [pos.isxcf] );
    lookup = trial.vs_filenames{ii};
    rows = find(cellfun(@(c) ischar(c) && strcmp(c, lookup), {pos.cmp}));
    lookup_T = trial.T_filenames{ii};
    rows_T = find(cellfun(@(c) ischar(c) && strcmp(c, lookup_T), {pos.cmp}));

    %if(sum(cell2mat({pos_red(rows).isxcf}))==16) %positions found
    %end
    trial.center_x{ii} = {pos(rows).center_x};
    trial.center_y{ii} = {pos(rows).center_y};
    trial.istarget{ii} = {pos(rows).istarget};
    trial.item{ii}     = {pos(rows).item};
%     pos_red = pos( [pos.isxcf] );
%     lookup = trial.vs_filenames{ii};
%     rows = find(cellfun(@(c) ischar(c) && strcmp(c, lookup), {pos_red.cmp}));
%     %if(sum(cell2mat({pos_red(rows).isxcf}))==16) %positions found
%     %end
%     trial.center_x{ii} = {pos_red(rows).center_x};
%     trial.center_y{ii} = {pos_red(rows).center_y};
%     trial.istarget{ii} = {pos_red(rows).istarget};
%     trial.item{ii}     = {pos_red(rows).item};

end

end

