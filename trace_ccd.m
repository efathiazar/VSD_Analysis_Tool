function out_trace=trace_ccd(ccd1,BW,do_normalize)
    %ccd1=handles.A{kk};
    if nargin < 2
        error('Not enough parameters')
    elseif nargin == 2
        do_normalize = false;
    end
    out_trace=zeros(1);
    if size(size(ccd1))==[1 3]
        [m n z]=size(ccd1);
        for i=1:z
            b(:,:)=ccd1(:,:,i);%+5);
            out_trace(i,1)=sum(b(BW==1))/size(find(BW==1),1);
        end
    elseif size(size(ccd1))==2
        out_trace(1,1)=sum(ccd1(BW==1))/size(find(BW==1),1);
    end
    if do_normalize
        if length(out_trace)>=20
            tt(:,1)=out_trace(10:20);
            means_s=sum(tt(:))/size(tt,1);
            out_trace(:,1)=out_trace(:,1)/means_s;
        end
    end
    %if ~isempty(s)
    %    handles.s{handles.Border_Num}=s;
    %end
    %plot(handles.axes2,s)
    
%     for i=1:size(dataw_a,3)
%         b(:,:)=dataw_a(:,:,i);%+5);
%         s(i,1)=sum(b(BW==1))/size(find(BW==1),1);
%     end

end