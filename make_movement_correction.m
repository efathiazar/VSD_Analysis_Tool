%%%%%%%%%%% This function gets one 3D data 
% load data and choose the trial you want to use
%%%%%%%%%%% and start and end frame as input and gives the original
%%%%%%%%%%% data and corrected version of data as output
%%% if ~start or ~end then start(s)=1, end=size(data,3)
%function [images,imagew,s,e]=make_movement_correction(image,reference_im,s,e)
function [images,imagew,s,e]=make_movement_correction(image,s,e)
%keyboard
%     for i1=1:nargin
%         if verLessThan('matlab', '8.3.1')
%         % -- Put code to run under MATLAB 8.3.0 and earlier here --
%             input_var.(genvarname(inputname(i1)))                = varargin{i1};
%         else
%         % -- Put code to run under MATLAB 8.3.1 and later here -- 
%             input_var.(matlab.lang.makeValidName(inputname(i1))) = varargin{i1};
%         end
%         input_var.varnames{i1}                               = inputname(i1);    
%     end
%     if ~isfield(input_var,'experiment')
%         warning('No experiment data available')
%     else
    if nargin==1
		s=1;
		e=size(image,3)-1;
		if e==0
			e=size(image,4);
		else
			e=size(image,3);
		end
	elseif nargin==2
		e=size(image,3)-1;
		if e==0
			e=size(image,4);
		else
			e=size(image,3);
        end
    end
	addpath('mex');
	%keyboard
	image=mat2gray(image);
%% correct movement from frame to frame and plot trace of the image

%     if nargin==1
        im11(:,:)=image(:,:,s+5); % 6th frame
%     elseif nargin==2
%         if size(reference_im,1) == size(image(:,:,s+5),1) &&...
%                 size(reference_im,2) == size(image(:,:,s+5),2)
%             im11 = reference_im;
%         else
%             fprintf('Ratio error of reference image')
%             im11(:,:)=image(:,:,s+5); % 6th frame
%         end
%     end
	for j=1:3
		im1(:,:,j)=im11;
	end
%%%%%%%%%%%%%%
	for i=s:e
		%if i==s;
			%im11(:,:)=image(:,:,s);
			%BW1=BW;
		im22(:,:)=image(:,:,i);
		%else
			%im11(:,:)=imagew(:,:,i+1-s);
			%BW1=BW;
		%	im22(:,:)=image(:,:,i+1);
		%end
		%im11=mat2gray(im11);
		%im22=mat2gray(im22);
		%[u,v,cert] = HierarchicalLK(im1, im2, 2, 1, 5, 1);
		for j=1:3
			%im1(:,:,j)=im11;
			im2(:,:,j)=im22;
		end
		warpI2=demflow(im1, im2);
		images(:,:,i+1-s)=rgb2gray(im2);
		imagew(:,:,i+1-s)=rgb2gray(warpI2);
	end
end


