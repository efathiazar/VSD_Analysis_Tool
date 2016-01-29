function vsd_trace_gui(varargin)
    %% function main_figure = hah_vsd_gui(experiment,boundary_tall,varargin)
    % obligatory input:
    %                   experiment -    struct with the 1dim-cell field  
    %                                   data_snapshots, where each cell is one
    %                                   2dim-Matrix (picture)
    %                   
    % output:
    %                   
    %
    if nargin >= 1
        for i1=1:nargin
            if verLessThan('matlab', '8.3.1')
            % -- Put code to run under MATLAB 8.3.0 and earlier here --
                input_var.(genvarname(inputname(i1)))                = varargin{i1};
            else
            % -- Put code to run under MATLAB 8.3.1 and later here -- 
                input_var.(matlab.lang.makeValidName(inputname(i1))) = varargin{i1};
            end
            input_var.varnames{i1}                               = inputname(i1);    
        end
        if ~isfield(input_var,'experiment')
            warning('No experiment data available')
        else
            plotchoice.medianC1 = 0;
            plotchoice.medianC2 = 0;
            plotchoice.medianC12 = 0;
            plotchoice.num_trials = 0;

            experiment = input_var.experiment;
            boundary_tall = resize_boundary( experiment );
            boundary_tall = boundary_tall(~cellfun(@isempty,boundary_tall));
            digital = input_var.experiment.analog_measurements{4,8};
            digital_time = 0;
            idx = 0;
            new = 1;
            current_BW = 1;
            for i = 1:length(digital)
                if digital(i)>4000 
                    if new == 1
                        idx = idx+1;
                        digital_time(idx,1) = i;
                    end
                    new = 0;
                else
                    new = 1;
                end
            end
%             if length(boundary_tall) < length(input_var.experiment.boundary)
%                 boundary_tall = input_var.experiment;
%             end
            if ~isfield(input_var,'main_figure')        
                input_var.main_figure = figure(1);
                hold off
            end
            set(input_var.main_figure, 'Units', 'normalized', 'Position', [0, 0, 1, 1]);
            colormap gray
            subplot(1,2,1)
            title('Snapshot with ROIs')
            if isfield(input_var.experiment,'data_snapshots')
                imagesc(input_var.experiment.data_snapshots{current_BW});
            else
                error('No Snapshots available')
            end
            axis image        
            if iscell(boundary_tall)
                [hp_bound,ht_bound] = plot_boundary( boundary_tall );
            end
            set(gcf, 'Units', 'normalized', 'Position', [0, 0, 1, 1]);
            hpopup_snap = uicontrol('Style', 'popup',...
                   'String', {1:length(input_var.experiment.data_snapshots)},...
                   'Units','normalized','Position',[.13 .1 .05 .02],...
                   'Callback', @setsnap);
            htext2 = uicontrol('Style','text',...
                    'String','Select Snapshot:',...
                    'Units','normalized','Position',[.13 .13 .06 .02]);
            hmain_figure = input_var.main_figure;
            %%        
            if ~isfield(input_var.experiment,'s_debleached')
                warning('No debleached data available')
            else
                
                
                set(hp_bound(current_BW),'Color','blue','LineWidth',3)
               hsub1 = subplot(2,2,2);
                plot(digital_time(5:end),input_var.experiment.traces_norm(5:end,4,current_BW)-1)
                hold on
                %hp1 = plot(digital_time(5:end),input_var.experiment.traces_norm_median(5:end,current_BW)-1);
              %  set(hp1,'Color','black','LineWidth',4)   
                %set(gca, 'YTickLabel', num2str(get(gca,'YTick')','%.4f'))  
                grid on
                yax1 = ylim;
             %   xaxis = xlim;
                title('Normalized traces of chosen ROIs')           
                hold off
               hsub2 = subplot(2,2,4);
                %keyboard
               plot(digital_time(5:end),input_var.experiment.s_debleached_one(:,4,current_BW))
               hold on
              % hp2 = plot(digital_time(5:end),input_var.experiment.s_debleached_one_median(:,current_BW));
                %set(hp2,'Color','black','LineWidth',4)
                %set(gca, 'YTickLabel', num2str(get(gca,'YTick')','%.4f'))
                hp2 = plot(digital_time(5:end),input_var.experiment.s_debleached_median_stim_cell1(:,current_BW),'-');
                set(hp2,'Color','black','LineWidth',2)
                hp3 = plot(digital_time(5:end),input_var.experiment.s_debleached_median_stim_cell2(:,current_BW),'-.');
                set(hp3,'Color','blue','LineWidth',2)
                hp4 = plot(digital_time(5:end),input_var.experiment.s_debleached_median_stim_cell12(:,current_BW),'--');
                set(hp4,'Color','red','LineWidth',2)
                grid on
                yax2 = ylim;
                legend([hp2,hp3,hp4],'Stim Cell 1','Stim Cell 2','Stim Both')
                title('Debleached traces of chosen ROIs')
               hold off


                 yaxis_new(1) = min([yax1(1) yax2(1)]);
                 yaxis_new(2) = max([yax1(2) yax2(2)]);  
%                 linkaxes([hsub1 hsub2])
%                 ylim(yaxis_new);
%                 %         axis([xaxis yaxis_new]);
%keyboard
        subplot(hsub1)
        % ylim([-0.03,0.03]);
         axis([2000 13000 yaxis_new]);
         subplot(hsub2)
         axis([2000 13000 yaxis_new]);
               
%                text('String','Chose ROI:','Position',[1000 150 90 40]);
                hpopup_debleach = uicontrol('Style', 'popup',...
                   'String', {1:length(input_var.experiment.BWs)},...
                   'Units','normalized','Position',[.5 .1 .05 .02],...
                   'Callback', @setdebleach);
               % % create a static text to show log messages
                htext = uicontrol('Style','text',...
                    'String','BW: 1','Units','normalized','Position',[.5 .07 .06 .02]);
                htext1 = uicontrol('Style','text',...
                    'String','Choose BW:',...
                    'Units','normalized','Position',[.5 .13 .06 .02]);

            end
        end
    end
    %keyboard    
    %% subfunction
    function setsnap(source,callbackdata)
        if verLessThan('matlab', '8.3.1')
        % -- Put code to run under MATLAB 8.3.0 and earlier here --
            val = get(source,'Value');
        else
        % -- Put code to run under MATLAB 8.3.1 and later here -- 
            val = source.Value;
        end
        subplot(1,2,1)
        title('Snapshot with ROIs')
        imagesc(input_var.experiment.data_snapshots{val});
        axis image
        if iscell(boundary_tall)
            [hp_bound,ht_bound] = plot_boundary( boundary_tall );
        end
    end    
    function setdebleach(source,callbackdata)
        if verLessThan('matlab', '8.3.1')
        % -- Put code to run under MATLAB 8.3.0 and earlier here --
            val = get(source,'Value');
        else
        % -- Put code to run under MATLAB 8.3.1 and later here -- 
            val = source.Value;
        end
        set(hp_bound(current_BW),'Color','red','LineWidth',2)
        current_BW = val;
        set(hp_bound(val),'Color','blue','LineWidth',3)
        set(htext, 'String', {['BW: ' num2str(val)]});
        subplot(hsub1);
        plot(digital_time(5:end),input_var.experiment.traces_norm(5:end,:,val)-1)
        hold on
        %hp1 = plot(digital_time(5:end),input_var.experiment.traces_norm_median(5:end,val)-1);
       % set(hp1,'Color','black','LineWidth',4)
        yax1 = ylim;
        %set(gca, 'YTickLabel', num2str(get(gca,'YTick')','%.4f'))
        grid on
        title('Normalized traces of chosen ROIs')
        hold off
        subplot(hsub2);
        plot(digital_time(5:end),input_var.experiment.s_debleached_one(:,:,val))
        hold on
        hp2 = plot(digital_time(5:end),input_var.experiment.s_debleached_median_stim_cell1(:,val),'-');
        set(hp2,'Color','black','LineWidth',2)
        hp3 = plot(digital_time(5:end),input_var.experiment.s_debleached_median_stim_cell2(:,val),'.-');
        set(hp3,'Color','blue','LineWidth',2)
        hp4 = plot(digital_time(5:end),input_var.experiment.s_debleached_median_stim_cell12(:,val),'--');
        set(hp4,'Color','red','LineWidth',2)
        yax2 = ylim;
        %set(gca, 'YTickLabel', num2str(get(gca,'YTick')','%.4f'))
        grid on        
        %-|-.|--|:
        legend([hp2,hp3,hp4],'Stim Cell 1','Stim Cell 2','Stim Both')
        title('Debleached traces of chosen ROIs')
         hold off
%   
          yaxis_new(1) = min([yax1(1) yax2(1)]);
          yaxis_new(2) = max([yax1(2) yax2(2)]);  
%         linkaxes([hsub1 hsub2])
        subplot(hsub1)
        % ylim([-0.03,0.03]);
         axis([2000 13000 yaxis_new]);
         subplot(hsub2)
         axis([2000 13000 yaxis_new]);
    end    
end