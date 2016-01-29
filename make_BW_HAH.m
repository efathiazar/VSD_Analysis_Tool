% boundary handling % let's you to draw a boundary on the image and
% displays the drawn boundary
function [BW ,position]=make_BW_HAH(image,h)
    if nargin>=1
        im=[];
        if size(size(image),2)==2
            im=image;
        elseif size(size(image),2)==3
            im(:,:)=image(:,:,1);
        else
            if size(size(image),2)==4
               im(:,:)=image(:,:,1,1); 
            end
        end
        imagesc(im),colormap gray
        axis square
    end
    if nargin<2
       h = imellipse;
    %    h = impoly;
    end
    position = wait(h);
    BW = createMask(h);
    [B,L] = bwboundaries(BW,'noholes');
    boundary = B{1};  
    hold on
    plot(boundary(:,2),boundary(:,1),'g','LineWidth',2)
    hold off
    delete(h) 
end