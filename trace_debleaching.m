function s_debleached = trace_debleaching(s1,s2)
    %% the main part of debleaching function
    % Variables:
    %    Input: s1 any of those traces
    %           s2 reference trace without stimulus
    %    Output: s_debleached the debleaching line
    %     p2 = polyfit(1:size(ss),ss(:,1)',1);
    %     f2 = polyval(p2,1:size(ss));
    %     %keyboard
    %     ss_debleached = ss-f2';
    
    tt(:,1)=s1(10:20);
    means1=sum(tt(:))/size(tt,1);
    s_norm(:,1)=s1(:,1)/means1;
    
%     popup_choose = uicontrol('Style', 'popup',...
%                    'String', {'Trial:', 1:size(input_var.experiment.analog_measurements,1)},...
%                    'Position', [0.6 0.15 0.05 0.03],...
%                    'Callback', @choose_debleach);
    
%     function choose_debleach(source,callbackdata)
%         if verLessThan('matlab', '8.3.1')
%         % -- Put code to run under MATLAB 8.3.0 and earlier here --
%             val = get(source,'Value');
%         else
%         % -- Put code to run under MATLAB 8.3.1 and later here -- 
%             val = source.Value;
%         end
%         h2 = subplot(1,2,2)
%         s2 =         

        tt2(:,1)=s2(10:20);
        means2=sum(tt2(:))/size(tt2,1);
        s2_norm(:,1)=s2(:,1)/means2;

        p1 = polyfit(1:size(s_norm)-4,s_norm(5:end,1)',1);
        f1 = polyval(p1,1:size(s_norm)-4);

        p2 = polyfit(1:size(s2_norm)-4,s2_norm(5:end,1)',1);
        f2 = polyval(p2,1:size(s2_norm)-4);

        s_debleached = s_norm(5:end)-f2';
%    end
end