function [v,a]=load_experiment(path,file_name)
    %% FUNCTION: load_experiment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % input:    path        - path where the data is stored
    %           file_name   - name of the currently processed experiment
    %           %% not here: Memory_Leak_Protection -  false/true
    % output:   a - experiment data (Trial-number,VSD-data,Analog-Data,Rest of v)
    %           v - all data of the last trial
    %           file(a,v)
    %
    % Authors: 	Elham Fathiazar
    % Edited:   Helge Ahrens
    % Date: 	02.03.2015
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % loads a vscope experiment and extracts the 
    switch nargin
        case 0
            [FileName,PathName]=uigetfile({'*.mat;*.xml'});
        case 1
            [FileName,PathName]=uigetfile(strcat(path,'\','*.xml')); 
        otherwise
            FileName=file_name;
            PathName=path;
    end

    % addpath for vscope_load
    addpath(PathName);
    %addpath(genpath('JEFG'));
    %addpath(genpath('mex'));

    %keyboard
if findstr( '.xml',FileName)
    %[str,remain]=strtok(FileName,'.');
    listing=dir(strcat(PathName,'0*.xml') );
    N=1;
    for i=1:size(listing,1)
        nam=listing(i).name;
        if findstr( '-',nam)
        else
            [str,~]=strtok(nam,'.');
            v=vscope_load({PathName, str2num(str)});
            v_last = v;
               %keyboard
                a{N,1}=str;
                if isfield(v,'ccd')
                    a{N,2}=v.ccd.dat;
                    v.ccd.dat = 'a{:,2}';
                end
                if isfield(v,'analog')
                    a{N,3}=v.analog.dat; % new line
                    v.analog.dat = 'a{:,3}';
                end
                a{N,4} = v;
%                if ~Memory_Leak_Protection
                    N=N+1;
%                end
    %% due to memory leak the ammount of loaded trials is to be restricted
    %             if Memory_Leak_Protection
    %                 filename = [path_new,'a',int2str(ii),'.mat'];
    %                 save(filename,'a','v');
    %                 clear a
    %                 N=1;
    %                 ii=ii+1;
    %             end
    %%
            end
           % keyboard
        end
    %    keyboard
    %     if ~Memory_Leak_Protection
    %         filename = [path_new,'a_full','.mat'];
    %         save(filename,'a','v','-v7.3');
    %     end
    else
        mode = struct('WindowStyle','non-modal',... 
                            'Interpreter','tex');
        h = errordlg('No *.xlm file found in specified folder.',...
            'Error', mode);
    end
   % keyboard
    v = v_last;
end