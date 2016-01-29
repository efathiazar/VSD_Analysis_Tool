%%%%%%%%%%%%%% apply movement correction on the result of load experiment
function [datas,dataw]=experiment_movement_correction(data,s,e);
addpath('D:\Analysis+Track\mex');
%load('a18c.mat');
datas{1}=[];
dataw{1}=[];
for i=1:size(data,1)
    if size(size(data{i,2}),2)>2
 %%%%%%%%%%%%%%% it is not optimum because it stes s and e everytime      
if nargin==1
    s=1;
    e=size(data{i,2},4)-1;
else if nargin==2
        e=size(data{i,2},4)-1;
    end
end
 %%%%%%%%%%%%%%%%%%%%%%     
 clear image images imagew
image(:,:,:)=data{i,2}(:,:,1,:);
[images,imagew]=make_movement_correction(image,s,e);
datas{i}=images;
dataw{i}=imagew;
    end
end
