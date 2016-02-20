function vsd_gui(varargin)
    %% function vsd_gui(experiment,data_norm,dataw)
    % obligatory input:
    %                   experiment -    struct with the 1dim-cell field  
    %                                   data_snapshots, where each cell is one
    %                                   2dim-Matrix (picture)
    %                   data_norm
    %                   dataw
    % output:
    %                   -
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
            choose_frame = 6;
            %boundary_tall = cell(1);
            experiment = input_var.experiment;
            boundary_tall = resize_boundary( experiment );
            %keyboard
            if ~iscell(boundary_tall)
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
                boundary_tall = cell(0);
            else
                boundary_tall = boundary_tall(~cellfun(@isempty,boundary_tall));
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            end
%             if length(boundary_tall) < length(input_var.experiment.boundary)
%                 boundary_tall = input_var.experiment.boundary;
%             end
            if ~isfield(input_var,'main_figure')        
                input_var.main_figure = figure(1);
                hold off
            end
            set(input_var.main_figure, 'Units', 'normalized', 'Position', [0, 0, 1, 1]);
            colormap gray
            set(gcf, 'Units', 'normalized', 'Position', [0, 0, 1, 1]);
            subplot(1,2,1)

            if isfield(input_var.experiment,'data_snapshots')
                imagesc(input_var.experiment.data_snapshots{1});
            else
                error('No Snapshots available')
            end
            axis image        
            if iscell(boundary_tall)
                [h.p_bound,h.t_bound] = plot_boundary( boundary_tall );
            end
            subplot(1,2,2)
            if isfield(input_var,'dataw')
                imagesc(input_var.dataw(:,:,choose_frame,1));
                [h.p_bound_small,h.t_bound_small] = plot_boundary( input_var.experiment.boundary );
            else
                error('No Raw video data available')
            end            
            axis square
            h.text1 = uicontrol('Style','text',...
                   'String','Select Snapshot:',...
                   'Units','normalized','Position',[.13 .13 .06 .02]);
            h.drop_snap = uicontrol('Style', 'popup',...
                   'String', {1:length(input_var.experiment.data_snapshots)},...
                   'Units','normalized','Position',[.13 .1 .05 .02],...
                   'Callback', @setsnap);
            h.text2 = uicontrol('Style','text',...
                   'String','Select Trial:',...
                   'Units','normalized','Position',[.5 .13 .06 .02]);
            h.drop_raw = uicontrol('Style', 'popup',...
                   'String', {1:size(input_var.data_norm,4)},...
                   'Units','normalized','Position',[.5 .1 .05 .02],...
                   'Callback', @setraw);
            main_figure = input_var.main_figure;
            %%
            subplot(1,2,2)
            if ~isfield(input_var,'data_norm') && ~isfield(input_var,'dataw')
                warning('No experiment data available')
            else
                
                h.button_BW = uicontrol('Style', 'pushbutton',...
                   'String', 'New BW',...
                   'Units','normalized','Position',[.5 .16 .05 .02],...
                   'FontSize',10,'FontWeight','bold',...
                   'BackgroundColor', [.9 .9 .9],...
                   'Callback', @addBW);
               
%                 if ~iscell(experiment.BWs)
%                     runtimevar.x = inputdlg('Enter number of new ROIS:','ROIS');
%                     runtimevar.n = round_HAH(str2num(runtimevar.x{:}),1);
%                     [experiment.BWs,experiment.positions,experiment.boundary]...
%                         = create_and_plot_BWs(input_var.data_norm,input_var.dataw,runtimevar.n,...
%                         experiment.BWs,experiment.boundary,size(experiment.data_snapshots{1}));
% %                     [experiment.BWs,experiment.positions,experiment.boundary,experiment.traces_norm,experiment.traces] = create_and_plot_BWs...
% %                         (input_var.data_norm,input_var.dataw,runtimevar.n,experiment.BWs,experiment.boundary,size(experiment.data_snapshots{1}));
%                 else
%                     runtimevar.BWs_tmp = experiment.BWs;
%                     runtimevar.boundary_tmp = experiment.boundary;
%                     runtimevar.x = inputdlg('Enter number of new ROIS:','ROIS');
%                     runtimevar.n = round_HAH(str2num(runtimevar.x{:}),1);
%                     [experiment.BWs,experiment.positions,experiment.boundary]...
%                         = create_and_plot_BWs(input_var.data_norm,input_var.dataw,runtimevar.n,...
%                         experiment.BWs,experiment.boundary,size(experiment.data_snapshots{1}));
% %                     [experiment.BWs,experiment.positions,experiment.boundary,experiment.traces_norm,experiment.traces] = create_and_plot_BWs...
% %                         (input_var.data_norm,input_var.dataw,runtimevar.n,experiment.BWs,experiment.boundary,size(experiment.data_snapshots{1}));
%                     clear('runtimevar.BWs_tmp');
%                 end
                h.button_trace= uicontrol('Style', 'pushbutton',...
                   'String', 'Trace All',...
                   'Units','normalized','Position',[.8 .1 .05 .02],...
                   'FontSize',10,'FontWeight','bold',...
                   'BackgroundColor', [.9 .9 .9],...
                   'Callback', @doTraceAll);
                h.button_save = uicontrol('Style', 'pushbutton',...
                   'String', 'Save Data',...
                   'Units','normalized','Position',[.9 .1 .05 .02],...
                   'FontSize',10,'FontWeight','bold',...
                   'BackgroundColor', [.9 .9 .9],...
                   'Callback', @save_Data);
%                 date_string = clock;
%                 save(['K:\Studium\VSD\p_t_cell\results\130410.new.',...
%                     num2str(date_string(1)),num2str(date_string(2)),...
%                     num2str(date_string(3)),'.analog_and_vsd.mat'],...
%                     'experiment','miss','idx','data_norm','dataw','analog')
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
        imagesc(input_var.experiment.data_snapshots{val});
        axis image
        if iscell(boundary_tall)
            [h.p_bound,h.t_bound] = plot_boundary( boundary_tall );
        end
    end
    function setraw(source,callbackdata)
        if verLessThan('matlab', '8.3.1')
        % -- Put code to run under MATLAB 8.3.0 and earlier here --
            val = get(source,'Value');
        else
        % -- Put code to run under MATLAB 8.3.1 and later here -- 
            val = source.Value;
        end
        subplot(1,2,2)
        imagesc(input_var.dataw(:,:,choose_frame,val));
        axis square
        if iscell(input_var.experiment.boundary)
            [h.p_bound_small,h.t_bound_small] = plot_boundary( input_var.experiment.boundary );
        end
    end

    function addBW(source,callbackdata)
        if verLessThan('matlab', '8.3.1')
        % -- Put code to run under MATLAB 8.3.0 and earlier here --
            val = get(source,'Value');
        else
        % -- Put code to run under MATLAB 8.3.1 and later here -- 
            val = source.Value;
        end
        if ~iscell(experiment.BWs)
           [experiment.BWs,experiment.positions,experiment.boundary]...
                = create_and_plot_BWs(input_var.data_norm,input_var.dataw,1,...
                experiment.BWs,experiment.boundary,size(experiment.data_snapshots{1}));
%                     [experiment.BWs,experiment.positions,experiment.boundary,experiment.traces_norm,experiment.traces] = create_and_plot_BWs...
%                         (input_var.data_norm,input_var.dataw,runtimevar.n,experiment.BWs,experiment.boundary,size(experiment.data_snapshots{1}));
        else
            runtimevar.BWs_tmp = experiment.BWs;
            runtimevar.boundary_tmp = experiment.boundary;
            [experiment.BWs,experiment.positions,experiment.boundary]...
                = create_and_plot_BWs(input_var.data_norm,input_var.dataw,1,...
                experiment.BWs,experiment.boundary,size(experiment.data_snapshots{1}));
%                     [experiment.BWs,experiment.positions,experiment.boundary,experiment.traces_norm,experiment.traces] = create_and_plot_BWs...
%                         (input_var.data_norm,input_var.dataw,runtimevar.n,experiment.BWs,experiment.boundary,size(experiment.data_snapshots{1}));
        end
    end

    function save_Data(source,callbackdata)
        if verLessThan('matlab', '8.3.1')
        % -- Put code to run under MATLAB 8.3.0 and earlier here --
            val = get(source,'Value');
        else
        % -- Put code to run under MATLAB 8.3.1 and later here -- 
            val = source.Value;
        end
%         prompt = 'Enter a keyword for the file to be saved:';
%         dlg_title = 'Filename Additional field';
%         num_lines = 1;
%         def = {'new'};
        experiment;
        data_norm = input_var.data_norm;
        dataw = input_var.dataw;
        %filename_prefix=inputdlg(prompt,dlg_title,num_lines,def);
         date_string = clock;
         uisave({'experiment','data_norm','dataw'},['K:\Studium\VSD\p_t_cell\results\130410.',...
                    num2str(date_string(1)),num2str(date_string(2)),...
                    num2str(date_string(3)),'.analog_and_vsd.mat']);
%                 save(['K:\Studium\VSD\p_t_cell\results\130410.',filename_prefix,'.',...
%                     num2str(date_string(1)),num2str(date_string(2)),...
%                     num2str(date_string(3)),'.analog_and_vsd.mat'],...
%                     'experiment','data_norm','dataw')
        msgbox('Data saved!')
    end

    function doTraceAll(source,callbackdata)
        if verLessThan('matlab', '8.3.1')
        % -- Put code to run under MATLAB 8.3.0 and earlier here --
            val = get(source,'Value');
        else
        % -- Put code to run under MATLAB 8.3.1 and later here -- 
            val = source.Value;
        end
        experiment = trace_all( experiment,input_var.data_norm,input_var.dataw );
    end

    function [BWs,positions,boundary,trial,trialw] =...
        create_and_plot_BWs(data_norm,dataw,n,BWs,boundary,snap_ratio)
        %% plot more cells, name them with 1, 2, 3 and show the number beside them
        %
        %
        if nargin==2
            n=1; % number of cells you want to show
        end
        trial=[];
        trialw=[];
        idxBWs = 0;
        BW2 = cell(size(BWs));
        if nargin>=4
            if ~iscell(BWs)
                BWs=cell(1,n);
            else
                idxBWs = size(BWs,2);
            end
            if nargin>=5
                if ~iscell(boundary)
                    boundary=cell(1,n);
                else
                    idxBWs = size(boundary,2);
                end
            else
                boundary=cell(1,n);
            end
            if nargin>=5
                if ndims(snap_ratio)== 1
                    snap_ratio = [snap_ratio snap_ratio];
                elseif ndims(snap_ratio) >= 2
                    snap_ratio = snap_ratio(1:2);
                end
            end
        else
            BWs=cell(1,n);
        end
        image=data_norm(:,:,6,1);
        num=size(data_norm,4);
        freq_vsd=100;
        t_vsd=5/freq_vsd:1/freq_vsd:(size(data_norm,3)-1)/freq_vsd;

        %% preallocation of matrices and vectors needed
        %% Data to be saved (maybe outside)
        for i4=idxBWs+1:idxBWs+n
            subplot(1,2,2)
            hold off
            %[BWs{i4},positions{i4}]=make_BW_HAH(image);
            [BWs{i4},positions{i4}]=make_BW_HAH();
            hold off
            [B,L] = bwboundaries(BWs{i4},'noholes');
            boundary{i4}=B{1};
            BW2{i4} = imresize(BWs{i4},snap_ratio);
            [B,L] = bwboundaries(BW2{i4},'noholes');
            boundary_tall{i4}=B{1};
            colors=['r' 'g'];
            if ~isempty(boundary_tall{i4})
                subplot(1,2,1)
                hold on
                plot_boundary( boundary_tall(i4),i4);
            end
        end
%         subplot(1,2,2)
%         hold on
%         for i4=1:idxBWs+n
%             if ~isempty(boundary{i4})
%                 plot(boundary{i4}(:,2),boundary{i4}(:,1),'g','LineWidth',2)
%                 if ~isempty(positions{i4})
%                     text(positions{i4}(1,1),positions{i4}(1,2),['\fontsize{14}\color{red} ' int2str(i4)],...
%                         'LineWidth',2,'HorizontalAlignment','left')
%                 end
%             drawnow
%             end
%         end
%         hold off
    end
end