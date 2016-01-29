% boundary handling % let's you to draw a boundary on the image and
% displays the drawn boundary
function [BW ,position,boundary]=edit_BW_HAH(image_handle,boundary)
    if nargin==2
        h = impoly(image_handle,[boundary(:,2),boundary(:,1)]);
        position = wait(h);
        BW = createMask(h);
        [B,L] = bwboundaries(BW,'noholes');
        boundary = B{1};  
        hold on
        plot(boundary(:,2),boundary(:,1),'g','LineWidth',2)
        hold off
        delete(h)
    else
        warning('Not enough parameters for editing the BW! Please provide image handle and boundary.')
    end    
end