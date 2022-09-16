function plot_scanpath(ET,trial,Exp,session_path,trials_to_plot)
%plot_scanpath
%
%% Example: Scanpath for an exemplary trial.
fdfactor = 0.05; % To adjust the size of the circles

screenXpixels = Exp.screenXpixels; %S10, S11
screenYpixels = Exp.screenYpixels; %S10, S11
subject       = Exp.subjname;

vs_filenames  = trial.vs_filenames;
T_filenames   = trial.T_filenames;
D1_filenames  = trial.D1_filenames;
D2_filenames  = trial.D2_filenames;
D3_filenames  = trial.D3_filenames;
D4_filenames  = trial.D4_filenames;

eyedata       = ET.VS.eyedata;

procdir       = Exp.procdir;

%trials_to_plot=[95 96 97 98 99]
for ii=1:length(trials_to_plot)
    tr = trials_to_plot(ii);

    Tpres         = trial.stim_present{5}(tr);
    respcorr      = trial.respcorr(tr);

    foo_target_count=0;
    for ii=1:210
        foo_target_count=foo_target_count+size(trial.istarget{ii},1);
    end

    labelpres     = numel(trial.item(tr)); %label with position exists?
    if Tpres==1 && labelpres==1
        str = trial.T_filenames(tr);
        Expr= '/';
        StartIndex = regexp(str,Expr);
        str_cell = char(str);
        str_cell = str_cell(cell2mat(StartIndex)+1:end);
        % D1
        lookup = str_cell;
        rows = find(cellfun(@(c) ischar(c) && strcmp(c, lookup), trial.item{tr}));
        index_target1_16 = rows;
        target_center_x = cell2mat(trial.center_x{tr}(index_target1_16));
        target_center_y = cell2mat(trial.center_y{tr}(index_target1_16));
        %change coordinate reference
        %keyboard
        %target_center_x = (Exp.screenXpixels - target_center_x)/2;
        %target_center_y = (Exp.screenYpixels - target_center_y)/2;
    else
        target_center_x = [];
        target_center_y = [];
    end




    figure; clf
    set(gcf,'Position',[Exp.screenXpixels/2 Exp.screenYpixels/2 Exp.screenXpixels/2 Exp.screenYpixels/2])
    set(gcf,'Color','w')
    set(gca,'XTick',[],'YTick',[]);
    set(gca,'visible','off')
    axes('Position',[0.05 0.26 0.55 0.55])
    set(gca,'XTick',[],'YTick',[]);
    set(gca,'visible','off')

    hold on
    % Image
    %keyboard
    A = imread(fullfile(session_path.images_folder,vs_filenames{tr}));
    %%
    %figure(7)
    xim=Exp.screenXpixels/2-size(A,2)/2;
    yim=Exp.screenYpixels/2-size(A,1)/2;
    imagesc(xim,yim,A)
    set(gca,'XTick',[],'YTick',[]);

    if respcorr ==1
        feedback_color = [0 1 0];
    else
        feedback_color = [1 0 0];
    end

    %set size of circle around each stimulus
    r = sqrt(2*40*40);
    t = linspace(0, 2*pi);

    % plot distractor positions as well
    Ndist=numel(trial.item{1});
    for iid=1:Ndist
        target_center_x_d = (Exp.screenXpixels - size(A,2))/2 + cell2mat(trial.center_x{tr}(iid));
        target_center_y_d = (Exp.screenYpixels - size(A,1))/2 + cell2mat(trial.center_y{tr}(iid));
        X = r*cos(t) + target_center_x_d;
        Y = r*sin(t) + target_center_y_d;
        patch(X,Y ,'k','FaceAlpha', 0.1, 'EdgeColor','none')
    end

    if length(target_center_x)==1
        %change coordinate reference frame
        target_center_x = (Exp.screenXpixels - size(A,2))/2 + target_center_x;
        target_center_y = (Exp.screenYpixels - size(A,1))/2 + target_center_y;
        X = r*cos(t) + target_center_x;
        Y = r*sin(t) + target_center_y;
        patch(X,Y ,feedback_color,'FaceAlpha', 0.25, 'EdgeColor',feedback_color)
    end

    %fprintf(' trial: %d\n',tr)
    %dbstop in plot_scanpath.m at 108 if tr==96
    if eyedata(tr).Nfix > 0
        % Scanpath
        plot(   eyedata(tr).samples(:,2), ...
            eyedata(tr).samples(:,3),'-','Color',[0 0 0.2],'LineWidth',2)

        % Fixations
        fixdur  = eyedata(tr).fixs(1:end,3)*fdfactor;
        xfix    = eyedata(tr).fixs(1:end,4);
        yfix    = eyedata(tr).fixs(1:end,5);
        nfix    = 1:length(fixdur);
        col = rainbow_colors(length(nfix));

        for i = 1:length(nfix)
            plot(xfix(i),yfix(i),'.','Color',col(i,:),'MarkerSize',ceil(fixdur(i)))
        end
        hold off
        set(gca,'XLim',[1 Exp.screenXpixels],'YLim',[1 Exp.screenYpixels])
        set(gca,'YDir','reverse','Visible','off')
        axes('Position',[0.05 0.82 0.55 0.01])
        hold on
        for i=1:length(nfix)
            hf=fill([i i+1 i+1 i i]-0.5,[1 1 0 0 1],col(i,:)); set(hf,'EdgeColor',col(i,:));
        end
        hold off
        if length(nfix)>4
            set(gca,'XLim',[0.5 length(nfix)-0.5],'XTick',1:5:length(nfix),'XAxisLocation','top')
        else
            set(gca,'XLim',[0.5 length(nfix)],'XTick',1:1:length(nfix),'XAxisLocation','top')
        end
        set(gca,'YLim',[0 1],'YTick',[])
        xlabel('Fixation Rank', 'fontsize',10,'fontweight','b','color','black');
        set(gca,                'fontsize',10,'fontweight','b')
        axes('Position',[0.10 0.08 0.50 0.15])
        hold on
        hx=plot(eyedata(tr).samples(:,1)/1000, eyedata(tr).samples(:,2),'-','Color',[.0 .0 .0],'LineWidth',2);
        hy=plot(eyedata(tr).samples(:,1)/1000, eyedata(tr).samples(:,3),'-','Color',[.5 .5 .5],'LineWidth',2);
        hold off
        xlim(eyedata(tr).samples([1 end],1)'/1000)
        ylim([0 max([Exp.screenXpixels,Exp.screenYpixels])])
        ylabel('eye position',  'fontsize',10,'fontweight','b','color','black');
        xlabel('time (secs)',   'fontsize',10,'fontweight','b','color','black');
        set(gca,                'fontsize',10,'fontweight','b')
        box on
        %%
        %         legend([hx hy],{'horizontal','vertical'},'Location','NorthEast')
        %D1
        t_img = imread(fullfile(session_path.images_folder,D1_filenames{tr}));
        ax2 = axes('Position',[0.05 0.92 0.075 0.075]);
        %ax2 = axes('Position',[0.6 0.8 0.075 0.075]);
        imagesc(t_img);
        set(gca,'XTick',[],'YTick',[]);
        %D2
        t_img = imread(fullfile(session_path.images_folder,D2_filenames{tr}));
        ax2 = axes('Position',[0.15 0.92 0.075 0.075]);
        imagesc(t_img);
        set(gca,'XTick',[],'YTick',[]);
        %D3
        t_img = imread(fullfile(session_path.images_folder,D3_filenames{tr}));
        ax2 = axes('Position',[0.25 0.92 0.075 0.075]);
        imagesc(t_img);
        set(gca,'XTick',[],'YTick',[]);
        %foo=sprintf('Scanpath S:%s, Trial: %d, \n Tpres: %d, Correct: %d',subject,tr, Tpres, respcorr );
        %title(foo);
        %D4
        t_img = imread(fullfile(session_path.images_folder,D4_filenames{tr}));
        ax2 = axes('Position',[0.35 0.92 0.075 0.075]);
        imagesc(t_img);
        set(gca,'XTick',[],'YTick',[]);
        %target
        t_img = imread(fullfile(session_path.images_folder,T_filenames{tr}));
        %ax2 = axes('Position',[0.7 0.50 0.15 0.15]);
        ax2 = axes('Position',[0.45 0.92 0.075 0.075]);
        img_target=imagesc(t_img);
        set(gca,'XTick',[],'YTick',[]);
    end

    annotation('rectangle',[0.445 0.915 0.085 0.085],'Color',feedback_color)
    foo={sprintf('Scanpath'),sprintf('S:%s, Trial:%d',subject,tr)};
    annotation('textbox',[0.60 0.9 0.075 0.075],'String',foo,'FitBoxToText','on');
    set(gca,'XTick',[],'YTick',[]);

    %
    if length(target_center_x)==1
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        outfile=fullfile(Exp.procdir,sprintf('scanpath_S%s_trial%d_fname%s',Exp.subjname,tr,vs_filenames{tr}));
        print(outfile,'-dpng','-r0')
    end
end
end

