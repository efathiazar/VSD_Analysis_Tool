function [v,a]=load_experiment_orig(path)
if nargin<1
[FileName,PathName]=uigetfile({'*.mat;*.xml' });
else
 [FileName,PathName]=uigetfile(strcat(path,'\','*.xml'));   
end
addpath(PathName);
% addpath for vscope_load
% addpath('D:\VSD_vscope\MATLAB\VSD_Experiments');
if findstr( '.xml',FileName)
    [str,remain]=strtok(FileName,'.');
    listing=dir(strcat(PathName,'0*.xml') );
    N=1;
    for i=1:size(listing,1)
        nam=listing(i).name;
        if findstr( '-',nam)
        else
            [str,remain]=strtok(nam,'.');
            v=vscope_load({PathName, str2num(str)});
           % if size(size(v.ccd.dat),2)>3
                a{N,1}=str;
                a{N,2}=v.ccd.dat;
                if isfield(v, 'analog')
                    a{N,3}=v.analog.dat; % new line
               % else
                %    a{N,3}=[]; % new line
                  
                end
                N=N+1;
            %end
        end
    end
end

