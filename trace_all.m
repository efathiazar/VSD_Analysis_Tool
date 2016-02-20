function [ experiment ] = trace_all( experiment,data_norm,dataw )
%% TRACE_ALL this function is to calculate the traces of all trials at ones
% input:
%       experiment -    BWs as a Cell field
%       dataw      -    movement corrected VSD data in 4dim-matrix
% output:
%       experiment -    traces       - traces without debleaching and
%                                      with movement correction
%                       traces_norm  - traces without debleaching and
%                                      with movement correction and
%                                      normalization
%                       s_debleached - traces with debleaching
%                                      with movement correction and
%                                      normalization
% ATTENTION - NOTE: if your input and output experiment variable is the same,
%                   these values will be overwritten and no backup is done
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msgbox('Trace All may take some time!')
    num.trials=size(dataw,4);
    num.frames = size(dataw,3);
    num.frames_debleached = num.frames-4;
    num.BWs=length(experiment.BWs);
    
    experiment.traces = zeros(num.frames,num.trials,num.BWs);
    experiment.traces_norm = experiment.traces;
    experiment.traces_no_movecorr = experiment.traces;
    experiment.s_debleached = zeros(num.frames_debleached,num.trials,num.BWs);
    experiment.s_debleached_one = experiment.s_debleached;
    experiment.s_debleached_mean = zeros(num.frames_debleached,num.BWs);
    experiment.s_debleached_one_mean = zeros(num.frames_debleached,num.BWs);
    experiment.traces_norm_mean = zeros(num.frames,num.BWs);
    
    for i4=1:num.BWs
            for j=1:num.trials
                    % tracing of non movement corrected without normalization           
                    dat(:,:,:)=data_norm(:,:,:,j);
                    experiment.traces_no_movecorr(:,j,i4)=...
                        trace_ccd(dat,experiment.BWs{i4},false);% mean
                    clear('dat')
                    % tracing of movement corrected without normalization           
                    datw(:,:,:)=dataw(:,:,:,j);
                    experiment.traces(:,j,i4)=...
                        trace_ccd(datw,experiment.BWs{i4},false);% mean

                    % tracing of movement corrected with normalization
                    experiment.traces_norm(:,j,i4)=...
                        trace_ccd(datw,experiment.BWs{i4},true);% mean
                    clear('datw')
            end
            experiment.traces_norm_mean(:,i4)=mean(experiment.traces_norm(:,:,i4),2);% mean
            experiment.traces_norm_median(:,i4)=median(experiment.traces_norm(:,:,i4),2);% median
    end

    num.analog = size(experiment.analog_measurements,1);
    analog.stim_sums = zeros(num.analog,2);
    analog.stim_median = zeros(num.analog,2);
    analog.stim_max = zeros(num.analog,2);
    analog.stim_cell = experiment.analog_measurements(:,3:4);   
    for ii = 1:num.analog
        analog.stim_max(ii,1) = max(analog.stim_cell{ii,1});
        analog.stim_max(ii,2) = max(analog.stim_cell{ii,2});
    end
    % Trials with stimulus
    analog.logic_idx = analog.stim_max>1;    
    % Trials with stimulus on both cells
    analog.idx_stim_cell12 = (sum(analog.logic_idx,2)==2);
    % Trials with stimulus on one cell
    analog.idx_stim_cell1 = ((analog.logic_idx(:,1)-analog.idx_stim_cell12)==1);
    analog.idx_stim_cell2 = ((analog.logic_idx(:,2)-analog.idx_stim_cell12)==1);
    
    % These are the indeces for trials without stimulus
    analog.logic_idx = analog.stim_max<1;
    % both cells no stim
    analog.logic_idx = (sum(analog.logic_idx,2)==2);
    % indeces of both cells no stim
    analog.logic_idx = find(analog.logic_idx==1);
    
    if numel(analog.logic_idx)==0        
        warning('There are no trials without stimulus!')
    end

    j2 = 1;
    for i4=1:num.BWs
            if numel(analog.logic_idx)>=j2
                for j=1:num.analog                
                    if j > analog.logic_idx(j2) && numel(analog.logic_idx)>j2
                        j2 = j2+1;
                    end
                    experiment.s_debleached(:,j,i4) =...
                        trace_debleaching(experiment.traces_norm(:,j,i4),...
                        experiment.traces_norm(:,analog.logic_idx(j2),i4));
                    experiment.s_debleached_one(:,j,i4) =...
                        trace_debleaching(experiment.traces_norm(:,j,i4),...
                        experiment.traces_norm(:,analog.logic_idx(1),i4));                    
                end
            else
                experiment.s_debleached(:,j,i4) = experiment.traces_norm(5:end,j,i4);
                experiment.s_debleached_warning =...
                    'There were no trials without stimulus for debleching!';
                for j=1:num.analog 
                    experiment.s_debleached_one(:,j,i4) =...
                            trace_debleaching(experiment.traces_norm(:,j,i4),...
                            experiment.traces_norm(:,1,i4));
                end
            end        
            experiment.s_debleached_mean(:,i4) = mean(experiment.s_debleached(:,:,i4),2);
            experiment.s_debleached_one_mean(:,i4) = mean(experiment.s_debleached_one(:,:,i4),2);
            experiment.s_debleached_median(:,i4) = median(experiment.s_debleached(:,:,i4),2);
            experiment.s_debleached_median_stim_cell1(:,i4) = median(experiment.s_debleached(:,analog.idx_stim_cell1,i4),2);
            experiment.s_debleached_median_stim_cell2(:,i4) = median(experiment.s_debleached(:,analog.idx_stim_cell2,i4),2);
            experiment.s_debleached_median_stim_cell12(:,i4) = median(experiment.s_debleached(:,analog.idx_stim_cell12,i4),2);
            experiment.s_debleached_one_median(:,i4) = median(experiment.s_debleached_one(:,:,i4),2);

            experiment.s_debleached_one_median_stim_cell1(:,i4) = median(experiment.s_debleached_one(:,analog.idx_stim_cell1,i4),2);
            experiment.s_debleached_one_median_stim_cell2(:,i4) = median(experiment.s_debleached_one(:,analog.idx_stim_cell2,i4),2);
            experiment.s_debleached_one_median_stim_cell12(:,i4) = median(experiment.s_debleached_one(:,analog.idx_stim_cell12,i4),2);
    end
    % experiment.traces = traces;
    % experiment.traces_norm = traces_norm;
    %end

    % 
    % figure_children = get(gcf,'Children')
    % get(figure_children,'Type')
    % set(figure_children(1),'Position',[1000 150 100 50])
    msgbox('Trace All done!')
end