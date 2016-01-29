function undo_plot(plot_handle)
%function undo_plot(h, n)
% UNDO_PLOT Undo last 'n' plotting operations.
%     if nargin == 1 
%       n = 1;
%     end
%     if n < 1
%       error('Cannot undo < 1 plotting operation!');
%     end
%     figure(h);
%     children = get(gca, 'children');
%     for i = 1 : n
%      delete(children(i));
      delete(plot_handle);
%    end
end

