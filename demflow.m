function warpI2=demflow(im2, im1)
%% im1 im2 to register first to be converted to another odor
% ebteda im1 be im2 register mishavad ta beshavad BW ra be digari tabdil
% kard
%addpath('mex');
alpha = 0.012;
ratio = 0.75;
%ratio = 0.9;
minWidth = 20;
nOuterFPIterations = 1;
nInnerFPIterations = 7;
nSORIterations =30;

para = [alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations];


% this is the core part of calling the mexed dll file for computing optical flow
% it also returns the time that is needed for two-frame estimation
% Coarse2FineTwoFrames: function to compute dense optical flow field in a 
%   coarse to fine manner
%tic;
% evalc is suppressing standard output of this function as long as it is
% not an error
[T,vx,vy,warpI2] = evalc('Coarse2FineTwoFrames(im2,im1,para);');
%[vx,vy,warpI2] = Coarse2FineTwoFrames(im2,im1,para);
%toc
%BW1=mat2gray(BW1);%
%[warpBW1,I]=warpFL(BW1,vx,vy);
%figure;imshow(im1);figure;imshow(warpI2);