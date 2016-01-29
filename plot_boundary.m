function [p_bound, t_bound] = plot_boundary( boundary,number )
%PLOT_BOUNDARY plots boundaries in active plot/subplot with numbering them
% p_bound = plot_boundary( boundary,number )
% input: boundary - cell with elements (n,2)
%        number - if only one boundary given, which one is it?
% ouput: handles of the plots
hold on
colors=['r' 'g'];
    for i=1:size(boundary,2)
        if ~isempty(boundary{i})
            p_bound(i) = plot(boundary{i}(:,2),boundary{i}(:,1),colors(1),'LineWidth',2);
            %randomize text position for better visibility
            rndRow = ceil(length(boundary{i})/(mod(rand*i,7)+1));
            col = boundary{i}(rndRow,2); 
            row = boundary{i}(rndRow,1);
            if size(boundary,2)>1
                number = i;
            end
            t_bound(i) = text(col+1, row-1, num2str(number));
            set(t_bound(i),'Color',colors(2),...
                'FontSize',14,'FontWeight','bold');
        end
    end
    drawnow
hold off
end

