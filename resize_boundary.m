function [ varargout ] = resize_boundary( varargin )
    %% RESIZE_BOUNDARY resiezes boundaries on the basis of given BWs and a ratio
    % obligatory input:
    %                   experiment -    struct with: 
    %                                   1dim-cell field data_snapshots, where each cell is one
    %                                   2dim-Matrix (picture)
    %                                   1dim-cell field BWs
    %keyboard
    if nargin >= 1
        %% initalize variables
        no_ratio = false;
        boundary = 0;
        idxBound = 0;
        %%%%%%%%%%%%%%%%%%%%%%
        for i2=1:nargin
            if verLessThan('matlab', '8.3.1')
            % -- Put code to run under MATLAB 8.3.0 and earlier here --
                input_var.(genvarname(inputname(i2)))               = varargin{i2};
            else
            % -- Put code to run under MATLAB 8.3.1 and later here -- 
                input_var.(matlab.lang.makeValidName(inputname(i2))) = varargin{i2};
            end            
            input_var.varnames{i2}                               = inputname(i2);    
        end
       %keyboard
        if ~isfield(input_var,'experiment')
            warning('No experiment data available')
            varargout = {false};
        else
            if isfield(input_var,'snap_ratio')
                if ndims(input_var.snap_ratio)==1
                    input_var.snap_ratio = [input_var.snap_ratio input_var.snap_ratio];
                elseif ndims(input_var.snap_ratio) >= 2
                    input_var.snap_ratio = input_var.snap_ratio(1:2);
                else
                    no_ratio = true;
                end
            elseif isfield(input_var,'experiment')
                if isfield(input_var.experiment,'data_snapshots')
                    if iscell(input_var.experiment.data_snapshots) && ...
                            size(input_var.experiment.data_snapshots,1)>=1
                        input_var.snap_ratio = size(input_var.experiment.data_snapshots{1});
                        input_var.snap_ratio = input_var.snap_ratio(1:2);
                    else
                        no_ratio = true;
                    end
                else
                    no_ratio = true;
                end
            else
                no_ratio = true;
            end
            if no_ratio% default ratio 512x512
                input_var.snap_ratio = [512,512];
            end
            if isfield(input_var.experiment,'BWs')
                if size(input_var.experiment.BWs,2) > 0
                    if ~isfield(input_var,'idxBWs')
                        input_var.idxBWs = size(input_var.experiment.BWs,2);
                    else
                        if input_var.idxBWs <= 0
                            input_var.idxBWs = 1;
                        end
                    end
                    boundary = cell(1,input_var.idxBWs);
                    %keyboard
                    for i3=1:input_var.idxBWs
                        if ~isempty(input_var.experiment.BWs{i3})
                            idxBound = idxBound+1;
                            BW2{i3} = imresize(input_var.experiment.BWs{i3},input_var.snap_ratio);
                            [B,L] = bwboundaries(BW2{i3},'noholes');
                            if ~isempty(B)
                                boundary{idxBound}=B{1};
                            else
                                idxBound = idxBound-1;
                            end
                        end
                    end
                end
                varargout = { boundary };
            else
                warning('No BWs available')
                %keyboard
                varargout = {false};
            end
        end
    end
end

