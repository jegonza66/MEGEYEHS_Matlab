function eyedata = fun_parse_eyedata(all,t_bgn,t_end)
% modified from fun_parse_eyedata
% adapted for right eye only
% read timing from Psychopy output ???
% MJI, 18.02.2020
% Expanded again to right and left eye
% May 20, 2022


msg     = [all.msg];
msgtime = [all.msgtime];
if strcmp(all.ojo,'LEFT')
    efix  = all.lefix;
    esac = all.lesac;
    ebli = all.lebli;
elseif strcmp(all.ojo,'RIGHT')
    efix  = all.refix;
    esac = all.resac;
    ebli = all.rebli;
else
    error('check which EYE was used during tracking')
end
    %t_bgn = msgtime(cellfun(@(x) ~isempty(strfind(x,'Task')),msg));
%t_end = msgtime(cellfun(@(x) ~isempty(strfind(x,'Retain')),msg));

eyedata = [];

fprintf('We are using fixations that could be partially included in the interval of interest\n')
fprintf('\ti.e. ending after the begin or starting before the end.\n')
for tr = 1:length(t_bgn)    
    % fix: t_start t_end dur x_medio y_medio pupila_medio
    indfix = (efix(:,2)>t_bgn(tr) & efix(:,1)<t_end(tr));
    Nfix = sum(indfix);
    if (Nfix>0)
        fixations   = nan(Nfix, size(efix,2));
        valid_fix   = find(indfix==1);
        for i = 1:length(valid_fix)
            fixations(i,:)  = efix(valid_fix(i),:);
        end
        
        indsac = (esac(:,2)>t_bgn(tr) & esac(:,1)<t_end(tr));
        Nsac = sum(indsac);
        saccades   = nan(Nfix, size(esac,2));
        valid_sac   = find(indsac==1);
        for i = 1:length(valid_sac)
            saccades(i,:)  = esac(valid_sac(i),:);
        end
        
        xpos = fixations(:,4);
        ypos = fixations(:,5);
            
        eyedata(tr).Nfix      = Nfix;
        eyedata(tr).fixs      = fixations;
        eyedata(tr).Nsac      = Nsac;
        eyedata(tr).sacs      = saccades;
        ind_trial             = all.samples(:,1)>t_bgn(tr) & all.samples(:,1)<t_end(tr);
        if(sum(ind_trial)==0)
            eyedata(tr).samples = [];
            eyedata(tr).fixs      = [];
            eyedata(tr).sacs      = [];
        else 
            eyedata(tr).samples   = all.samples(ind_trial,:);
            eyedata(tr).samples(:,1) = eyedata(tr).samples(:,1)-eyedata(tr).samples(1,1);
        end
    else
        indsac = (esac(:,2)>t_bgn(tr) & esac(:,1)<t_end(tr));
        Nsac = sum(indsac);   
        eyedata(tr).Nfix      = Nfix;
        eyedata(tr).fixs      = [];
        eyedata(tr).Nsac      = Nsac;
        eyedata(tr).sacs      = [];
        eyedata(tr).samples   = [];
    end
end
end