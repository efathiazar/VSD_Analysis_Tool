function VSD_analysis(varargin)
    %% VSD Load and proceed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This script is extracting the 4D-Matrix of VSD experiments and stores
    % them in a folder called result in one upper hirachy than the source
    % experiment data
    %
    % Authors: 	Helge Ahrens
    % Date: 	02.03.2015
    % Edited:   10.03.2015 % added memory leak switch variable
    %           11.02.2016 % added movementcorrection frame 6 (Oliver Kuehn)
    %           11.02.2016 % changed save(folders.FileName,'a','v','-v7.3');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin >= 1
        for i1=1:nargin
            if verLessThan('matlab', '8.3.1')
            % -- Put code to run under MATLAB 8.3.0 and earlier here --
                input_var.(genvarname(inputname(i1)))                = varargin{i1};
            else
            % -- Put code to run under MATLAB 8.3.1 and later here -- 
                input_var.(matlab.lang.makeValidName(inputname(i1))) = varargin{i1};
            end
            input_var.varnames{i1}                                   = inputname(i1);    
        end
    end
    %% Workspace initialization
    %    clear all
    %    close all
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% For debugging
    %      dbstop if error
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Bool-toggle Variables
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        do_bool.do_load_all_experiments             = false; 
        % if none of both, please choose experiment %
        % OR %%%
        do_bool.do_load_experiment_mat              = false;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        do_bool.do_save_experiment_mat              = true;
        do_bool.do_save_stimuli                     = true;
        % movement correction needs quite a long time
        do_bool.do_movement_correction_preprocess   = true;
        do_bool.do_movement_correction_on_reference = true;
        do_bool.do_movement_correction              = true;
        do_bool.do_create_new_ROIS                  = false;
        do_bool.do_load_existing_ROIS               = true;
        do_bool.do_load_existing_ROIS_plot          = false;

        %% not implemented yet %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        do_bool.do_create_plots                     = false;    %% maybe more precise? 
        do_bool.noise_test                          = false;    %% not implemented yet
        do_bool.Memory_Leak_Protection              = false;    %% if true then each trial 
                                                     %of an experiment is stored 
                                                     %as extra file in the results 
                                                     %folder (should be concidered 
                                                     %to be usedin case the machine has 
                                                     %less than 8GBphysical memory

        %keyboard
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% initial variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FOLDERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %folders.initial_experiments_folder = ['D:\VSD_vscope\MATLAB\VSD_Experiments'];
        %folders.initial_experiments_folder = ['K:\Studium\VSD'];
        folders.initial_experiments_folder =...
                ['D:\VSD_vscope\MATLAB\VSD_Experiments\p_t_cell'];
        if ~exist(folders.initial_experiments_folder,'dir')
            folders.initial_experiments_folder = ...
                    uigetdir(pwd,'Please choose a valid experiment folder');
    %keyboard
            if ~ischar(folders.initial_experiments_folder)
                error('No folder specified!')
            end

        end
        %folders.PathName = [folders.initial_experiments_folder,filesep,'2013'];
        folders.PathName = [folders.initial_experiments_folder];
        folders.PathName_results =...
            [folders.initial_experiments_folder,filesep,'results',filesep];
        folders.experiment_folders=dir(strcat(folders.PathName)); %,'0*.xml'));

        % VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        BWs                             = {};
        positions                       = {};
        boundary                        = {};
        traces_norm                     = {};
        traces                          = {};    
        data_snapshots                  = {};
        a                               = {};
        v                               = 0;
        analog_data                     = {};
        experiment.BWs                  = 0;
        experiment.positions            = 0;
        experiment.boundary             = 0;
        experiment.traces_norm          = 0;
        experiment.traces               = 0;
        experiment.data_snapshots       = 0;
        experiment.v_all                = 0;
        experiment.analog_data          = 0;
        miss.info                       =...
            ['In case there was a ratio mismatch between the trial data,',...
            ' these excluded mismatches are stored in this structure,',...
            ' there indeces can be found in the struct idx.ratiomismatch'];
        runtimevar.a                    = 0;
        runtimevar.v                    = 0;
        runtimevar.v_err                = 0;
        runtimevar.sleep_timer          =.01;
        runtimevar.VSD_size_over_200MB  = {}; %% in case of Memory_Leak_Protection to be used
    %% ROIS (regions of interest) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        runtimevar.n = 20; % number of ROIS (regions of interest)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% create subfolder for results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %path_idx = find(folders.PathName == filesep,2,'last');
        %path_new = [folders.PathName(1:path_idx(1)),'results'];
            %...,filesep,folders.PathName(path_idx(1)+1:path_idx(2))];
        if ~exist(folders.PathName_results,'dir')
            mkdir(folders.PathName_results);
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% adding a function recources to PATH variable
      %addpath(genpath('C:\Users\Helge\Documents\MATLAB\Experiment_VSD'));
      addpath(genpath('JEFG'));
      addpath(genpath('mex'));
      addpath(genpath('fit_ellipse'));
      addpath(genpath('piotr_toolbox'));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    %% Loading all experiments after eachother and saving the results
    if do_bool.do_load_experiment_mat
    %     if do_bool.do_load_experiment_mat_list
    %         folders.listing=dir(strcat(folders.nam,'*.mat') );
    %         if ~isempty(listing)
    %         end
    %     else
        [folders.FileName,folders.PathName]=uigetfile(strcat(folders.PathName,'\','*.mat'));    
        load([folders.PathName,folders.FileName]);
    %     end
    else
        if do_bool.do_load_all_experiments
            for i=1:size(folders.experiment_folders,1)
                folders.nam=[folders.PathName,filesep,folders.experiment_folders(i).name,filesep];
                %keyboard
                if ~strcmp(folders.nam,[folders.PathName,filesep,'routines',filesep]) &&...
                        ~strcmp(folders.nam,[folders.PathName,filesep,'results',filesep]) &&...
                        ~strcmp(folders.nam,[folders.PathName,filesep,'scripts',filesep]) &&...
                        ~strcmp(folders.nam,[folders.PathName,filesep,'_scripts',filesep]) &&...
                        ~strcmp(folders.nam,[folders.PathName,filesep,'_results',filesep]) &&...
                        ~strcmp(folders.nam,[folders.PathName,filesep,'matlab',filesep]) &&...
                        ~strcmp(folders.nam,[folders.PathName,filesep,'.',filesep]) &&...
                        ~strcmp(folders.nam,[folders.PathName,filesep,'..',filesep])
                    %% Load the experiment        
                    folders.listing=dir(strcat(folders.nam,'0*.xml') );
                    if ~isempty(folders.listing)
                        %[runtimevar.v,runtimevar.a] = load_experiment(folders.nam,listing(1).name,Memory_Leak_Protection);
                        [runtimevar.v,runtimevar.a] = load_experiment(folders.nam,folders.listing(1).name);
                        folders.FileName = [folders.PathName_results,folders.experiment_folders(i).name,'.mat'];
                        a = runtimevar.a;
                        v = runtimevar.v;
                        save(folders.FileName,'a','v','-v7.3');
                    end
                end
            end
            do_bool.do_load_experiments = false;
        else
            [runtimevar.v,runtimevar.a]=load_experiment_orig(folders.initial_experiments_folder);  
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% copy information from the old format to the new
    if iscell(BWs)
        experiment.BWs         = BWs;
        clear('BWs');
    end
    if iscell(positions)
        experiment.positions   = positions;
        clear('positions');
    end
    if iscell(boundary)
        experiment.boundary    = boundary;
        clear('boundary');
    end
    if ~isempty(traces_norm)
        experiment.traces_norm = traces_norm;
        clear('traces_norm');
    end
    if ~isempty(traces)
        experiment.traces      = traces;
        clear('traces');
    end
    if ~isempty(data_snapshots)
        experiment.data_snapshots = data_snapshots;
        clear('data_snapshots');
    end
    if ~isempty(analog_data)
        experiment.analog_data = analog_data;
        clear('analog_data');
    end
    if ~isempty(a)
        runtimevar.a = a;
        clear('a');
    end
    if isstruct(v)
        runtimevar.v = v;
        clear('v');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if iscell(experiment.v_all)&& ~isstruct(runtimevar.v)
        runtimevar.v = experiment.v_all{kk};
    end
    if isstruct(runtimevar.v)
        if isfield(runtimevar.v,'info')
    %         folders.experiment_file = [folders.PathName,filesep,'results',filesep,...
    %             runtimevar.v.info.expt,'_processed.mat'];
            folders.experiment_file = [folders.PathName,...
                runtimevar.v.info.expt,'_processed.mat'];
        else
    %         folders.experiment_file = [folders.PathName,filesep,'results',filesep,...
    %             'runtimevar.v.info.expt.nonexist','_processed.mat'];
            folders.experiment_file = [folders.PathName,...
                'runtimevar.v.info.expt.nonexist','_processed.mat'];
        end
    else
        error('Information of VSD project not available!')
    end

    %% Determine number of snapshots and real measurements %%%%%%%%%%%%%%%%%%%%
    if iscell(runtimevar.a)
        idx.num = 1;
        idx.snap = 1;
        for kk = 1:size(runtimevar.a,1)
          if ndims(runtimevar.a{kk,2}) == 4
            idx.measurement(idx.num) = kk;%excluding the snapshots 
            idx.analog(idx.num) = size(runtimevar.a{kk,3},2);
            idx.num = idx.num+1;
          else
              if ~isempty(runtimevar.a{kk,2})
                idx.snapshots(idx.snap) = kk;
                idx.snap = idx.snap+1;
              end
          end
        end
    %elseif 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Movement Correction % Preprocess %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if do_bool.do_movement_correction_preprocess
      %%%% normalization (all data together) % gave errors caused by snapshots
      %keyboard
      if iscell(runtimevar.a)
        %experiment.data_snapshots=zeros(size(runtimevar.a{idx.snapshots(1),2},1),...
        %    size(runtimevar.a{idx.snapshots(1),2},2),length(idx.snapshots));

        %keyboard
        %data_norm=mat2gray(runtimevar.v.ccd.dat(:,:,1,:));
        %dat=zeros(size(data_norm,1),size(data_norm,2),size(data_norm,4));
        data=zeros(size(runtimevar.a{idx.measurement(1),2},1),...
            size(runtimevar.a{idx.measurement(1),2},2),size(runtimevar.a{idx.measurement(1),2},4)...
            ,length(idx.measurement));
        experiment.analog_data = zeros(size(runtimevar.a{idx.measurement(1),3},1),...
            size(runtimevar.a{idx.measurement(1),3},2),length(idx.measurement));
        idx.miss=0;
        idx.miss_analog=0;
        idx.num=0;
        for i=idx.measurement
            if size(size(runtimevar.a{i,2}),2)>2
                idx.num=idx.num+1;
                if size(data(:,:,:,idx.num),1) == size(runtimevar.a{i,2}(:,:,1,:),1) &&...
                    size(data(:,:,:,idx.num),2) == size(runtimevar.a{i,2}(:,:,1,:),2) &&...
                    size(data(:,:,:,idx.num),3) == size(runtimevar.a{i,2}(:,:,1,:),4)

                    data(:,:,:,idx.num)=runtimevar.a{i,2}(:,:,1,:);
                else
                    idx.miss_data = idx.miss_data+1;                
                    miss.data_with_ratiomismatch{idx.miss_data} = runtimevar.a{i,2}(:,:,1,:);
                    idx.ratiomismatch.data(idx.miss_data) = i;
                end

                if size(experiment.analog_data(:,:,idx.num),1) == size(runtimevar.a{i,3}(:,:),1)&&...
                    size(experiment.analog_data(:,:,idx.num),2) == size(runtimevar.a{i,3}(:,:),2)

                    experiment.analog_data(:,:,idx.num) = runtimevar.a{i,3}(:,:);
                else
                    idx.miss_analog = idx.miss_analog+1;
                    miss.analog_with_ratiomismatch{idx.miss_analog} = runtimevar.a{i,3}(:,:);
                    idx.ratiomismatch.analog(idx.miss_analog) = i;
                end
            end
        end
        idx.snap=0;
        if ~iscell(experiment.data_snapshots)
            experiment.data_snapshots=cell(1,length(idx.snapshots));
            for i=idx.snapshots
                %experiment.data_snapshots(:,:,idx.snap) = runtimevar.a{i,2};
                idx.snap=idx.snap+1;
                experiment.data_snapshots{idx.snap} = runtimevar.a{i,2};
            end        
        end
        %    keyboard
        if size(runtimevar.a,2)==4
            experiment.v_all = runtimevar.a(:,4);
        end
        runtimevar.a = 0;
      %      keyboard
        data_norm = mat2gray(data);
        %keyboard
        clear data;    
        dataw=zeros(size(data_norm,1),size(data_norm,2),...
          size(data_norm,3),size(data_norm,4));
        if ~do_bool.do_movement_correction
            dataw=data_norm;
            experiment.dataw_info = 'no movement correction';
        end
      else
          warning('There is probably no movement correction possible, no raw VSD data available. (variable: runtimevar.a)')
      end   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    %% Movement Correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if do_bool.do_movement_correction
        % Olivers addition  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % makes movement correction for reference picture nr. 6 %%%%%%%%
         if do_bool.do_movement_correction_on_reference
             references(:,:,:) = data_norm(:,:,6,:);
             references2 = references;
             references2(:,:,6:size(references,3)+5) = references2;
             [ref_norm(:,:,:),ref_w(:,:,:),~,~] = make_movement_correction(references2);
             ref_w2 = ref_w(:,:,6:size(ref_norm,3));
             data_norm(:,:,6,:) = ref_w2;
             clear references references2 ref_norm ref_w2 ref_w
         end
        %keyboard
            %reference_im(:,:)=data_norm(:,:,6,1); % 6th frame of the first trial
         for i=1:size(data_norm,4)
            %clear images imagew
            %image(:,:,:)=data_norm(:,:,:,i);
            %[datas(:,:,:,i),dataw(:,:,:,i),~,~] = make_movement_correction(data_norm);
            [data_norm(:,:,:,i),dataw(:,:,:,i),~,~] = make_movement_correction(data_norm(:,:,:,i));
            %[data_norm(:,:,:,i),dataw(:,:,:,i),~,~] = make_movement_correction(data_norm(:,:,:,i),reference_im);
            %data_norm/datas is the same as original with normalization
            %dataw movement corrected data
         end
         dataw=gray2mat_ccd(dataw);
         data_norm=gray2mat_ccd(data_norm);
         experiment.dataw_info = 'with movement correction';
        %keyboard
    end
    if do_bool.do_load_existing_ROIS
    %% LOAD EXISTING ROIS % experimental %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %keyboard
        do_bool.do_overwrite_BWs = true;
        if (~iscell(experiment.BWs)&& ~iscell(experiment.boundary)) || do_bool.do_overwrite_BWs
            idx.idxBWs = 0;
            BWs_bkp = experiment.BWs;
            experiment.BWs = cell(size(runtimevar.v.rois));
            BW2 = cell(size(runtimevar.v.rois));
            experiment.boundary = cell(size(runtimevar.v.rois));
            for kk = 1:size(runtimevar.v.rois,2)
                if size(runtimevar.v.rois{kk},2)==2
                    idx.idxBWs = idx.idxBWs+1;
                    runtimevar.ellipse = fit_ellipse(runtimevar.v.rois{kk}(:,1),runtimevar.v.rois{kk}(:,2));
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % calculate boundaries and BWs for the fitted ellipse %%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    runtimevar.ellipse.cos_phi = cos( runtimevar.ellipse.phi );
                    runtimevar.ellipse.sin_phi = sin( runtimevar.ellipse.phi );
                    runtimevar.ellipse.R = ...
                        [ runtimevar.ellipse.cos_phi runtimevar.ellipse.sin_phi;...
                        -runtimevar.ellipse.sin_phi runtimevar.ellipse.cos_phi ];
                    runtimevar.ellipse.theta_r         = linspace(0,2*pi);
                    runtimevar.ellipse.ellipse_x_r     = ...
                        runtimevar.ellipse.X0 + runtimevar.ellipse.a*...
                        cos( runtimevar.ellipse.theta_r );
                    runtimevar.ellipse.ellipse_y_r     = ...
                        runtimevar.ellipse.Y0 + runtimevar.ellipse.b*...
                        sin( runtimevar.ellipse.theta_r );
                    experiment.boundary{idx.idxBWs} = ...
                        (runtimevar.ellipse.R * [runtimevar.ellipse.ellipse_x_r;...
                        runtimevar.ellipse.ellipse_y_r])';
                    BW2{idx.idxBWs} = poly2mask(experiment.boundary{idx.idxBWs}(:,1),experiment.boundary{idx.idxBWs}(:,2),...
                        size(experiment.data_snapshots{1,1},1),size(experiment.data_snapshots{1,1},2));
                    experiment.BWs{idx.idxBWs} = imresize(BW2{idx.idxBWs},[size(data_norm,1),size(data_norm,2)]);

                elseif size(runtimevar.v.rois{kk},2)==5
                    idx.idxBWs = idx.idxBWs+1;
                    %keyboard
                    experiment.boundary{idx.idxBWs} = fit_ellipse2(runtimevar.v.rois{kk});
                    BW2{idx.idxBWs} = poly2mask(experiment.boundary{idx.idxBWs}(:,1),experiment.boundary{idx.idxBWs}(:,2),...
                        size(experiment.data_snapshots{1,1},1),size(experiment.data_snapshots{1,1},2));
                    experiment.BWs{idx.idxBWs} = imresize(BW2{idx.idxBWs},[size(data_norm,1),size(data_norm,2)]);
                end
            end
            %keyboard
            experiment.boundary = experiment.boundary(~cellfun(@isempty,experiment.boundary));
        else
           idx.idxBWs = size(experiment.BWs,2);
           idx.idxboundary = size(experiment.boundary,2);
        end  
        if do_bool.do_load_existing_ROIS_plot
            imagesc(experiment.data_snapshots{1,1})
            colormap gray
            for kk = 1:size(experiment.boundary,2)
                if ~isempty(experiment.boundary{kk})
                    hold on
                    plot(experiment.boundary{kk}(:,1),experiment.boundary{kk}(:,2),'g','LineWidth',2)
                end
            end
            hold off
        end
       % keyboard
    end
    %keyboard
    if do_bool.do_create_new_ROIS
    %% CREATE ROIS % From here the tracing can be done %%%%%%%%%%%%%%%%%%%%%%%%
    %  Debleaching has to be reworked though
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % define the number of ROIS (regions of interest) in the VARIABLES section
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % trace incl. debleaching is included
    % vsd_gui(experiment,data_norm,dataw);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %keyboard
    end
    %% Save data in dedicated files, one for only the analog recording with 
    % experiment data and traces but without raw videa data
    if do_bool.do_save_stimuli
        if ~isstruct(idx)
            idx.measurement = 1:size(data_norm,4);
            idx.analog      = size(experiment.analog_data,2);        
        end
        experiment.analog_measurements  = cell(length(idx.measurement),max(idx.analog));
        experiment.info                 = cell(length(idx.measurement),1);    
        experiment.traces               = cell(length(idx.measurement),1);
        for kk = 1:length(idx.measurement)
            for ii = 1:max(idx.analog)
                if iscell(runtimevar.a)
                    if size(runtimevar.a{idx.measurement(kk),3},2)>=ii                
                        experiment.analog_measurements{kk,ii} = runtimevar.a{idx.measurement(kk),3}(:,ii);
                    end
                else
                    if size(experiment.analog_data,2)>=ii                
                        experiment.analog_measurements{kk,ii} = experiment.analog_data(:,ii,kk);
                    end                
                end
                if isfield(runtimevar.v,'analog')
                    experiment.info{kk}.analog_info       = runtimevar.v.analog.info;
                else
                    runtimevar.v_err = 1;
                end
                if isfield(runtimevar.v,'digital')
                    experiment.info{kk}.digital           = runtimevar.v.digital;
                else
                    runtimevar.v_err = 1;
                end
                if isfield(runtimevar.v,'ccd')
                    experiment.info{kk}.ccd_info          = runtimevar.v.ccd.info;
                else
                    runtimevar.v_err = 1;
                end
                if isfield(runtimevar.v,'rois')
                    experiment.info{kk}.rois              = runtimevar.v.rois;
                else
                    runtimevar.v_err = 1;
                end
                if isfield(runtimevar.v,'info')
                    experiment.info{kk}.info              = runtimevar.v.info;
                else
                    runtimevar.v_err = 1;
                end
                if isfield(runtimevar.v,'settings')
                    experiment.info{kk}.settings          = runtimevar.v.settings;
                else
                    runtimevar.v_err = 1;
                end
                if runtimevar.v_err == 1;
                    experiment.error_info{kk} =...
                        'Some Info data for the trial was not found';
                    runtimevar.v_err = 0;
                end
    %             experiment.BWs         = experiment.BWs;
    %             experiment.positions   = experiment.positions;
    %             experiment.boundary    = experiment.boundary;
    %             experiment.traces_norm = experiment.traces_norm;
    %             experiment.traces      = experiment.traces;
            end
        end    
        date_string = clock;    
        save([folders.experiment_file,'.',num2str(date_string(1)),num2str(date_string(2))...
            ,num2str(date_string(3)),'.analog_and_vsd.mat'],'experiment','miss',...
            'idx','data_norm','dataw')
    % In case the data_norm and dataw shall be later in the experiment
    % structure the data cann be manipulated like this:
    %     experiment = rmfield(experiment,'data_norm');
    %     experiment = rmfield(experiment,'dataw');
        save([folders.experiment_file,'.',num2str(date_string(1)),num2str(date_string(2))...
            ,num2str(date_string(3)),'.analog.mat'],'experiment','miss','idx')        
    end
    if do_bool.do_save_experiment_mat
        date_string = clock;
        save([folders.experiment_file,'.',num2str(date_string(1)),num2str(date_string(2))...
            ,num2str(date_string(3)),'.bkp.mat'])
    end
    %% end
    % dbquit
    fprintf('VSD Analysis function closed. Thank you for choosing flying with VSD airlines ;)\n')
end
