function s=trace(ccd1,BW)
%ccd1=handles.A{kk};
s=zeros(1);
if size(size(ccd1))==[1 3]
[m n z]=size(ccd1);
for i=1:z-20
    b(:,:)=ccd1(:,:,i+5);%+5);
    s(i,1)=sum(b(BW==1))/size(find(BW==1),1);
end
elseif size(size(ccd1))==2
    s(1,1)=sum(ccd1(BW==1))/size(find(BW==1),1);
end
%if ~isempty(s)
%    handles.s{handles.Border_Num}=s;
%end
%plot(handles.axes2,s)