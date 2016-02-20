gif_main_figure = figure(1);
set(gif_main_figure, 'Units', 'normalized', 'Position', [0, 0, 1, 1]);
for i = 1:1:size(data_norm,4)

     imagesc(data_norm(:,:,6,i));
     axis square 
     colormap gray
     % gif utilities
     drawnow;
     frame = getframe(1);
     im = frame2im(frame);
     [imind,cm] = rgb2ind(im,256);
     outfile = 'references.gif';

     % On the first loop, create the file. In subsequent loops, append.
     if i==1
         imwrite(imind,cm,outfile,'gif','DelayTime',.3,'loopcount',1);
     else
        imwrite(imind,cm,outfile,'gif','DelayTime',.3,'writemode','append');
     end

end
