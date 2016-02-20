image_sequence(:,:,:) = dataw(:,:,6,:);
for i = 1:1:size(image_sequence,3)
    colormap(gray)
    imagesc(image_sequence(:,:,i));
 
    % gif utilities
    drawnow;
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    outfile = 'references_corrected.gif';
 
    % On the first loop, create the file. In subsequent loops, append.
    if i==1
        imwrite(imind,cm,outfile,'gif','DelayTime',1,'loopcount',1);
    else
        imwrite(imind,cm,outfile,'gif','DelayTime',1,'writemode','append');
    end
 
end