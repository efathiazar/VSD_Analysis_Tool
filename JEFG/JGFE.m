function varargout = JGFE(varargin)
% JGFE M-file for JGFE.fig
%      JGFE, by itself, creates a new JGFE or raises the existing
%      singleton*.
%
%      H = JGFE returns the handle to a new JGFE or the handle to
%      the existing singleton*.
%
%      JGFE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in JGFE.M with the given input arguments.
%
%      JGFE('Property','Value',...) creates a new JGFE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before JGFE_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to JGFE_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help JGFE

% Last Modified by GUIDE v2.5 18-Sep-2014 18:44:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @JGFE_OpeningFcn, ...
                   'gui_OutputFcn',  @JGFE_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before JGFE is made visible.
function JGFE_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to JGFE (see VARARGIN)
%global  A Border_Num borders BWs v
%Dir=uigetdir;
handles.current_s=[];%dar border va newborder baraye zakhire kardane akharin trace rasm shode dar axes2
handles.BW=cell(1);
handles.Border_Num=0; % shows the Number of ploted borders to save accuordingly; it used in New_border
handles.borders=cell(1);% It is for saving borders; it is used in New_border
handles.trace_Num=0; % it saves the number of traces (s)
handles.trace=cell(1);% It is for saving traces of the mean values of pixels inside borders so we don't need to calculate them each time; 
handles.Analog={};
%it is used in New_border, trace_plot
path=get(handles.path,'String');
[FileName,PathName]=uigetfile({strcat(path,'*.mat;*.xml') });
if findstr( '.mat',FileName)
v=load(strcat(PathName,FileName));
%handles.v=load(FileName);
%v=load(FileName);
%a=fieldnames(handles.v);
a=fieldnames(v);
%dat1 = getfield(handles.v,a{1});
%handles.v=
handles.A=cell(1);
aa=a{1};
if strcmp(aa,'borders')
    handles.borders=v.borders;
    handles.BW=v.BWs;
    handles.trace=v.traces;
    handles.A=v.Trials(1,:);
    if size(v.Trials)>1
    b=v.Trials(2,:);
    set(handles.Trial_List,'String',[b']);
    end
    Border_Num=size(handles.borders,2);
    handles.Border_Num=Border_Num; 
    handles.trace_Num=size(handles.trace,2);
    set(handles.border,'String',[1:Border_Num]');
    handles.trace_Num=size(v.traces,2);
    set(handles.List_trace,'String',v.traces(1,:)');
else if findstr(aa,'dat')
        for i=1:size(a)
    %handles.v=getfield(handles.v,a{i});
    dat1=getfield(v,a{i});
    handles.A{i}=dat1.ccd.dat;  %this part is based on data structure
    set(handles.Trial_List,'String',[a]);
        end
    end
end
end
if findstr( '.xml',FileName)
    [str,remain]=strtok(FileName,'.');
    v=vscope_load({PathName, str2num(str)});
    if isfield(v,'analog')
    handles.Analog{1}=v.analog.dat;
    plot(handles.axes_Analog,handles.Analog{1}(:,2))
    else
        handles.Analog{1}=[];
    end
    
    set(handles.Analog_Channels,'String',size(handles.Analog{1},2))
  
    handles.A{1}(:,:,:)=v.ccd.dat(:,:,1,:);  %this part is based on data structure
    a{1,1}='dat1';
    set(handles.Trial_List,'String',[a]);      
end
%dat1=getfield(v,a{1});
%handles.A=mat2gray(dat1.ccd.dat);
%[m n]=size(a);
%for i=1:m
%    if ~isstruct(getfield(v,a{i}))     
%       v = rmfield(v, a{i});
%    end   
%end
set(handles.Analog_Num,'String',2)

set(handles.Farme_Num,'String',1)
set(handles.Trial,'String',get(handles.Trial_List,'Value'))
% v is the name of loaded experiment
%ccd1=dat1.ccd.dat;
%A=mat2gray(ccd1);


%handles.v=v;
axes(handles.axes1),imagesc(handles.A{1}(:,:,1)),colormap gray
%k=get(handles.Trial_List,'String');

% Choose default command line output for JGFE
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes JGFE wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Executes on button press in Load_Trial.
function Load_Trial_Callback(hObject, eventdata, handles)
% hObject    handle to Load_Trial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Dir=uigetdir;
path=get(handles.path,'String');
[FileName,PathName]=uigetfile({strcat(path,'*.mat;*.xml') });
if findstr( '.xml',FileName)
    [str,remain]=strtok(FileName,'.');
    v=vscope_load({PathName, str2num(str)});
    k=size(handles.A,2);
    handles.A{k+1}(:,:,:)=v.ccd.dat(:,:,1,:);%this part is based on data structure
    if isfield(v,'analog')
    handles.Analog{k+1}=v.analog.dat;
    else
        handles.Analog{k+1}=[];
       end
    b{1,1}=strcat('dat',num2str(k+1));
    a=get(handles.Trial_List,'String');
    set(handles.Trial_List,'String',[a;b]);
      
end
if findstr( '.mat',FileName)
v=load(strcat(PathName,FileName));
%handles.v=load(FileName);
%v=load(FileName);
%a=fieldnames(handles.v);
a=fieldnames(v);
%dat1 = getfield(handles.v,a{1});
%handles.v=
%handles.A=cell(1);
aa=a{1};
     if findstr(aa,'dat')
        for i=1:size(a)
    %handles.v=getfield(handles.v,a{i});
        dat1=getfield(v,a{i});
        k=size(handles.A,2);
    handles.A{k+1}=dat1.ccd.dat;  %this part is based on data structure
    if isfield(v,'analog')
    handles.Analog{k+1}=v.analog.dat;
    
    plot(handles.axes_Analog,handles.Analog{k+1}(:,2))
    else
        handles.Analog{k+1}=[];
    end
    set(handles.Analog_Channels,'String',size(handles.Analog{k+1},2))
     b{1,1}=strcat('dat',num2str(k+1));
    a=get(handles.Trial_List,'String');
   set(handles.Trial_List,'String',[a;b]);
        end
    end
end

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Farme_Num_Callback(hObject, eventdata, handles)
% hObject    handle to Farme_Num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global A
% Hints: get(hObject,'String') returns contents of Farme_Num as text
%        str2double(get(hObject,'String')) returns contents of Farme_Num as a double
Fr_Num = get(hObject,'String');
axes(handles.axes1)
Trial_Num=get(handles.Trial_List,'Value');
imagesc(handles.A{1,Trial_Num}(:,:,str2num(Fr_Num))),colormap gray
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Farme_Num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Farme_Num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in Next_Frame.
function Next_Frame_Callback(hObject, eventdata, handles)
% hObject    handle to Next_Frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global A
Fr_Num = get(handles.Farme_Num,'String');
set(handles.Farme_Num,'String',str2num(Fr_Num)+1)
axes(handles.axes1)
Trial_Num=get(handles.Trial_List,'Value');
imagesc(handles.A{1,Trial_Num}(:,:,str2num(Fr_Num)+1)),colormap gray
guidata(hObject,handles)

% --- Executes on button press in Previous_Frame.
function Previous_Frame_Callback(hObject, eventdata, handles)
% hObject    handle to Previous_Frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global A
Fr_Num = get(handles.Farme_Num,'String');
set(handles.Farme_Num,'String',str2num(Fr_Num)-1)
axes(handles.axes1)
Trial_Num=get(handles.Trial_List,'Value');
imagesc(handles.A{1,Trial_Num}(:,:,str2num(Fr_Num)-1)),colormap gray
guidata(hObject,handles)

% --- Executes on button press in Import_List.
function Import_List_Callback(hObject, eventdata, handles)
% hObject    handle to Import_List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName]=uigetfile('borders');
load(strcat(PathName,FileName)); 
handles.borders=borders;
handles.BW=BWs;
handles.Border_Num=size(borders,2);
%save borders borders
% neveshtane import dar halate adi va na ezafe kardan be list
set(handles.border,'String',[1:size(borders,2)]');
guidata(hObject,handles)

% --- Executes on button press in New_border.
function New_border_Callback(hObject, eventdata, handles)
% hObject    handle to New_border (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Border_Num A borders BWs ccd1 v
boundary=[];
%Tri = get(handles.Trial,'String');
%kk=str2double(Tri);
kk = get(handles.Trial_List,'Value');
Fr_Num = get(handles.Farme_Num,'String');
%axes(handles.axes1)
figure,imagesc(handles.A{kk}(:,:,str2num(Fr_Num))),colormap gray
str=get(handles.popupmenu1,'String');
val=get(handles.popupmenu1,'Value');
%str{val};
%switch str{val}
 %   case 'elliptic'
 %       s=1;   
  %  case 'plygonal'
   %     s=2;
%end
 if val==1       
        h = imellipse;
 else
     h=impoly;
 end
        position = wait(h);
        BW = createMask(h);
%edg=edge(BW,'canny');
%save BW BW
%save edg edg
%A=handles.A;
%save A A
%[fe1,fe2]=find(edg);
%boundary = bwtraceboundary(BW,[fe1(1), fe2(1)],'NE');
[B,L] = bwboundaries(BW,'noholes');
%imagesc(mat2gray(A(:,:,1))), colormap gray
%hold on
%for k = 1:length(B)
    boundary = B{1};
%    plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
%end
SA=get(handles.border,'String');
if iscell(SA)
%if SA(1)=='Listbox'
    handles.Border_Num=0;
else
    handles.Border_Num=size(SA,1);
end
% when there is new boundry in the following code it adds 1 to Border_Num
%updating borders, BWs and Border_Num
if ~isempty(boundary)
    handles.Border_Num=handles.Border_Num+1;
    handles.borders{handles.Border_Num}=boundary;
    handles.BW{handles.Border_Num}=BW;
    %borders=[borders;boundary];
%    BWs=[BWs;BW];
end
close Figure 1
Border_Num=handles.Border_Num;
Fr_Num = get(handles.Farme_Num,'String');
k=get(handles.border,'String');
set(handles.border,'String',[1:Border_Num]')
axes(handles.axes1),imagesc(handles.A{kk}(:,:,str2num(Fr_Num))),colormap gray
hold on
%save boundary boundary
plot(boundary(:,2),boundary(:,1),'g','LineWidth',2)
%text(position(1,1),position(1,2),['\fontsize{12}\color{green}\leftarrow ' int2str(Border_Num)],...
%    'HorizontalAlignment','left')
text(position(1,1),position(1,2),['\fontsize{12}\color{red}\leftarrow ' int2str(Border_Num)],...
    'HorizontalAlignment','left')
hold off
% adding to the border_List(Listbox)
ccd1=handles.A{kk};
[m n z]=size(ccd1);
s=trace(ccd1,BW);
%for i=1:z%-10
   % b(:,:)=ccd1(:,:,i);%+5);
  %  s(i,1)=sum(b(BW==1));%c
 % %  imh(i,1)=imhist(b(BW==1));
%end
% Naming the trace- trace is s, Name of trace is Name_s
st = get(handles.Trial_List,'String');

% dar vaghe dar inja bar asase trace ha zakhire shode va atelaate avali va
% akhari faghat baraye takmile etelaate trace ast. khob pas BW ro ham
% yakhire mikonam
if ~isempty(s)
    handles.trace_Num=handles.trace_Num+1;
    Name_s{1,1}=strcat(st{kk},'_Border','-',num2str(handles.Border_Num));% in ghalate, bayad ba tavajoh be shomareye border numgozari kone
    handles.trace{1,handles.trace_Num}=Name_s;% name trace inja zakhire mishavad
    handles.trace{2,handles.trace_Num}=s;% khode trace inja zakhire mishavad
    %handles.s{3,handles.Border_Num}=boundary;%khode boundary ro zakhire mikone, baraye ine ke az in trace be mahale border dastrasi peida konim vali baraye in nist ke etelaate digari az BW be dast biaiad. banabar in BW ra zakhire nemikonim.
    %handles.s{4,handles.Border_Num}=BW; 
    Contents=get(handles.List_trace,'String');
set(handles.List_trace,'String',[Contents;Name_s])
end

plot(handles.axes2,s)
handles.current_s=s;
%save s s
guidata(hObject,handles)

% --- Executes on button press in Export_List.
function Export_List_Callback(hObject, eventdata, handles)
% hObject    handle to Export_List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName]=uiputfile(strcat('Borders\',datestr(now,1),num2str(hour(now)),num2str(minute(now)),'.mat'));
borders=handles.borders;
BWs=handles.BW;
save(strcat(PathName,FileName),'borders','BWs'); 
guidata(hObject,handles)


% --- Executes on button press in delet_border.
function delet_border_Callback(hObject, eventdata, handles)
% hObject    handle to delet_border (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Border_Nu =get(handles.border,'Value');%-1;
Border_Numb =get(handles.border,'String');
%Border_Num=str2double(Border_Nu);
%Tri = get(handles.Trial,'String');
%kk=str2double(Tri);
kk=get(handles.Trial_List,'Value');
Fr_Num = get(handles.Farme_Num,'String');
%axes(handles.axes1)
%borders=handles.borders;
if Border_Nu(1,1)==1
    handles.borders={handles.borders{2:end}};
else if Border_Nu(1,1)==size(Border_Numb,2)
        handles.borders={handles.borders{1:Border_Numb(1,1)-1}};
    else 
   %     s=handles.borders;
   %     save s s
       handles.borders={handles.borders{1:Border_Nu(1,1)-1},handles.borders{Border_Nu(1,1)+1:end}};
    end
end
%Border_Numb=
%handles.border =set(handles.border,'String',);
%k=get(handles.border,'String');
%set(handles.border,'String',[k;Border_Num])
axes(handles.axes1)
imagesc(handles.A{kk}(:,:,str2num(Fr_Num))),colormap gray
if size(Border_Nu,2)>1
hold on
for i=2:size(Border_Nu,2)
%save boundary boundary
b=Border_Nu(1,i);
boundary=handles.borders{b-1};
plot(boundary(:,2),boundary(:,1),'g','LineWidth',2)
%text(position(1,1),position(1,2),['\fontsize{12}\color{green}\leftarrow ' int2str(Border_Num)],...
%    'HorizontalAlignment','left')
text(boundary(1,2),boundary(1,1),['\fontsize{12}\color{green}' int2str(b-1) '\rightarrow '],...
    'HorizontalAlignment','right')
end
hold off
axes(handles.axes2)
plot(handles.axes2,handles.s{Border_Nu(1,2)-1}) 
hold on
for i=2:size(Border_Nu,2)
%save boundary boundary
b=Border_Nu(1,i);
sk=handles.s{b-1};
plot(sk)
end
hold off
axes(handles.axes1)
end
handles.Border_Num=handles.Border_Num-1;
Border_Num=handles.Border_Num;
set(handles.border,'String',[1:Border_Num]','Value',1)
guidata(hObject,handles)

% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in border.
function border_Callback(hObject, eventdata, handles)
% hObject    handle to border (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%contents = get(hObject,'String');
%save contents
%adding new Member: con=[contents;'me'];
% Hints: contents = get(hObject,'String') returns border contents as cell array
%        contents{get(hObject,'Value')} returns selected item from border
Border_Num =get(hObject,'Value');%-1;
%Border_Num=str2double(Border_Nu);
%Tri = get(handles.Trial,'String');
%kk=str2double(Tri);
kk = get(handles.Trial_List,'Value');
Fr_Num = get(handles.Farme_Num,'String');
axes(handles.axes1)
imagesc(handles.A{kk}(:,:,str2num(Fr_Num))),colormap gray
hold on
for i=1:size(Border_Num,2)
%save boundary boundary
b=Border_Num(1,i);
%b=i;
boundary=handles.borders{1,b};
%save boundary boundary
plot(boundary(:,2),boundary(:,1),'g','LineWidth',2)
%text(position(1,1),position(1,2),['\fontsize{12}\color{green}\leftarrow ' int2str(Border_Num)],...
%    'HorizontalAlignment','left')
text(boundary(1,2),boundary(1,1),['\fontsize{12}\color{red}' int2str(b) '\rightarrow '],...
    'HorizontalAlignment','right')
end
hold off
axes(handles.axes2)
ccd1=handles.A{kk};
s=trace(ccd1,handles.BW{Border_Num(1,1)});
%plot(handles.axes2,handles.s{Border_Num(1,1)}) 
plot(handles.axes2,s)
handles.current_s=s;
%hold on
%for i=1:size(Border_Num,2)
%%save boundary boundary
%b=Border_Num(1,i);
%sk=trace(ccd1,handles.BW{b});
%%sk=handles.s{b};
%plot(sk)
%end
%%save sk sk
%handles.current_s=sk;
%hold off
axes(handles.axes1)

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function border_CreateFcn(hObject, eventdata, handles)
% hObject    handle to border (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Enter.
function Enter_Callback(hObject, eventdata, handles)
% hObject    handle to Enter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in update_List.
function update_List_Callback(hObject, eventdata, handles)
% hObject    handle to update_List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Load_Experiment.
function Load_Experiment_Callback(hObject, eventdata, handles)
% hObject    handle to Load_Experiment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global A v
% v ham global ast ta beshe A ro bahash taghyir dad masalan dar
% Trial_Number, inja az A global estefade nakardim vali chon taghyiresh
% dadim va momkene jahaye dge azash estefade konim
handles.A={};
handles.Analog={};
handles.current_s=[];%dar border va newborder baraye zakhire kardane akharin trace rasm shode dar axes2
handles.BW=cell(1);
handles.Border_Num=0; % shows the Number of ploted borders to save accuordingly; it used in New_border
handles.borders=cell(1);% It is for saving borders; it is used in New_border
handles.trace_Num=0; % it saves the number of traces (s)
handles.trace=cell(1);% It is for saving traces of the mean values of pixels inside borders so we don't need to calculate them each time; 
%it is used in New_border, trace_plot
path=get(handles.path,'String');
[FileName,PathName]=uigetfile({strcat(path,'*.mat;*.xml') });
if findstr( '.mat',FileName)
v=load(strcat(PathName,FileName));
%handles.v=load(FileName);
%v=load(FileName);
%a=fieldnames(handles.v);
a=fieldnames(v);
%dat1 = getfield(handles.v,a{1});
%handles.v=
handles.A=cell(1);
aa=a{1};
if strcmp(aa,'borders')
    handles.borders=v.borders;
    handles.BW=v.BWs;
    handles.trace=v.traces;
    handles.A=v.Trials(1,:);
    if size(v.Trials)>1
    b=v.Trials(2,:);
    set(handles.Trial_List,'String',[b'],'val',1);
    end
    Border_Num=size(handles.borders,2);
    handles.Border_Num=Border_Num; 
    handles.trace_Num=size(handles.trace,2);
    set(handles.border,'String',[1:Border_Num]','val',1);
    handles.trace_Num=size(v.traces,2);
    set(handles.List_trace,'String',v.traces(1,:)','val',1);
else if findstr(aa,'dat')
        for i=1:size(a)
    %handles.v=getfield(handles.v,a{i});
    dat1=getfield(v,a{i});
    handles.A{i}=dat1.ccd.dat;  %this part is based on data structure
    set(handles.Trial_List,'String',[a],'val',1);
        end
    end
end
end
if findstr( '.xml',FileName)
   % [str,remain]=strtok(FileName,'.');
   % v=vscope_load({PathName, str2num(str)});
   % handles.A{1}(:,:,:)=v.ccd.dat(:,:,1,:); %this part is based on data structure
   % a{1,1}='dat1';
   % set(handles.Trial_List,'String',[a],'val',1);
   
   
   addpath(PathName);
    [str,remain]=strtok(FileName,'.');
    listing=dir(strcat(PathName,'0*.xml') );
    N=1;
    for i=1:size(listing,1)
        nam=listing(i).name;
        if findstr( '-',nam)
        else
            [str,remain]=strtok(nam,'.');
            v=vscope_load({PathName, str2num(str)});         
            if size(size(v.ccd.dat),2)>3
                handles.A{N}(:,:,:)=v.ccd.dat(:,:,1,:);  %this part is based on data structure 
                if isfield(v,'analog')
                handles.Analog{N}=v.analog.dat;
                plot(handles.axes_Analog,handles.Analog{1}(:,2))
                else
                    handles.Analog{N}=[];
                end
                set(handles.Analog_Channels,'String',size(handles.Analog{N},2))
                
                a{N,1}=str2num(str);
               % a{N,1}=str;
                N=N+1;
            end
        end
    end
end   


%dat1=getfield(v,a{1});
%handles.A=mat2gray(dat1.ccd.dat);
%[m n]=size(a);
%for i=1:m
%    if ~isstruct(getfield(v,a{i}))     
%       v = rmfield(v, a{i});
%    end   
%end
set(handles.Trial_List,'String',[a]);
set(handles.Farme_Num,'String',1)
set(handles.Trial,'String',get(handles.Trial_List,'Value'),'val',1)
% v is the name of loaded experiment
%ccd1=dat1.ccd.dat;
%A=mat2gray(ccd1);
axes(handles.axes1)
%handles.v=v;
imagesc(handles.A{1}(:,:,1)),colormap gray

%k=get(handles.Trial_List,'String');

% Choose default command line output for JGFE
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes JGFE wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.



function Trial_Callback(hObject, eventdata, handles)
% hObject    handle to Trial_Number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global A v
% Hints: get(hObject,'String') returns contents of Farme_Num as text
%        str2double(get(hObject,'String')) returns contents of Farme_Num as a double
Tri = get(hObject,'String');
kk=str2double(Tri);
%save Trial % natije be shekle '2' az jemse ab ast
%kk=str2double(Trial);

set(handles.Trial_List,'Value',str2double(Tri));
%a=fieldnames(handles.v);
%save a a
%clear   A
%save kk kk
%da = getfield(handles.v,a{kk});
%save da da
%ccd1=da.ccd.dat;
%handles.A=da.ccd.dat;
%clear da
%save ccd1 ccd1
%g=mat2gray(ccd1);
%handles.A=g;
%save g g
%axes(handles.axes1)
Fr_Num=get(handles.Farme_Num,'String');
Analog_Num=get(handles.Analog_Num,'String');
if ~isempty(handles.Analog{kk})
plot(axes_Analog,handles.Analog{kk}(:,str2num(Analog_Num)))

end
set(handles.Analog_Channels,'String',size(handles.Analog{kk},2))
axes(handles.axes1)
imagesc(handles.A{kk}(:,:,str2num(Fr_Num))),colormap gray
% Hints: get(hObject,'String') returns contents of Trial_Number as text
%        str2double(get(hObject,'String')) returns contents of Trial_Number as a double
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Trial_Number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Trial_Number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Trial_List.
function Trial_List_Callback(hObject, eventdata, handles)
% hObject    handle to Trial_List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Trial_List contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Trial_List
kk= get(hObject,'Value');
set(handles.Trial,'string',num2str(kk));
Fr_Num=get(handles.Farme_Num,'String');
axes(handles.axes1)
set(handles.axes1, 'Xlim',[0 1],'Ylim',[0 1])
imagesc(handles.A{kk}(:,:,str2num(Fr_Num))),colormap gray
Analog_Num=get(handles.Analog_Num,'String');
if ~isempty(handles.Analog{kk})
   
plot(handles.axes_Analog,handles.Analog{kk}(:,str2num(Analog_Num)))
end
 set(handles.Analog_Channels,'String',size(handles.Analog{kk},2))
% Hints: get(hObject,'String') returns contents of Trial_Number as text
%        str2double(get(hObject,'String')) returns contents of Trial_Number as a double
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Trial_List_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Trial_List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in average.
function average_Callback(hObject, eventdata, handles)
% hObject    handle to average (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%v=handles.v;
%a=fieldnames(handles.v);
%[m n]=size(a);
[m n]=size(handles.A);
%k=getfield(handles.v,a{1});
%Ave=k.ccd.dat;
Ave=handles.A{1};
%for i=2:m
for i=2:n
    Ave=handles.A{i}+Ave;
end

%    dat1 = getfield(handles.v,a{i});
%    Ave=dat1.ccd.dat+Ave;
%end
%Ave=Ave./n;
handles.A{n+1}=Ave;
%  Trial=handles.Trial;
%v=set(v,)
%handles.v=setfield(handles.v,'Avr',Avr);
a=get(handles.Trial_List,'String');
set(handles.Trial_List,'String',[a;'Ave']);
b=[a;Ave];
%save Ave Ave
guidata(hObject,handles)
% --- Executes on button press in Show_all_borders.
function Show_all_borders_Callback(hObject, eventdata, handles)
% hObject    handle to Show_all_borders (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Border_N =get(hObject,'String');% it is a {} which the first one is NAN
%Border_Nu=Border_N{2:end};
%Border_Num=str2double(Border_Nu);
%Tri = get(handles.Trial,'String');
%kk=str2double(Tri);
%Fr_Num = get(handles.Farme_Num,'String');
%axes(handles.axes1)
%imagesc(handles.A{kk}(:,:,str2num(Fr_Num))),colormap gray
%hold on
%save boundary boundary
%boundary=handles.borders{Border_Num};
%plot(boundary(:,2),boundary(:,1),'g','LineWidth',2)
%text(position(1,1),position(1,2),['\fontsize{12}\color{green}\leftarrow ' int2str(Border_Num)],...
%    'HorizontalAlignment','left')
%text(boundary(1,2),boundary(1,1),['\fontsize{12}\color{green}' int2str(Border_Num) '\rightarrow '],...
%    'HorizontalAlignment','right')
%hold off
%plot(handles.axes2,handles.s{Border_Num})

% --- Executes on button press in Load_Border.
function Load_Border_Callback(hObject, eventdata, handles)
% hObject    handle to Load_Border (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName]=uigetfile('borders');
load(strcat(PathName,FileName));
Border_Num=handles.Border_Num;
for i=1:size(borders,2)
    handles.borders{Border_Num+i}=borders{i};
    handles.BW{Border_Num+i}=BWs{i};
    handles.Border_Num=handles.Border_Num+1;
end
    

Border_Num=handles.Border_Num;
%handles.Border_Num=size(borders,2);
%save borders borders
% neveshtane import dar halate adi va na ezafe kardan be list
set(handles.border,'String',[1:Border_Num]');
guidata(hObject,handles)

% --- Executes on button press in Hold_Borders.
function Hold_Borders_Callback(hObject, eventdata, handles)
% hObject    handle to Hold_Borders (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Hold_Borders



function order_Callback(hObject, eventdata, handles)
% hObject    handle to order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of order as text
%        str2double(get(hObject,'String')) returns contents of order as a double


% --- Executes during object creation, after setting all properties.
function order_CreateFcn(hObject, eventdata, handles)
% hObject    handle to order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in debleach_polyfit.

function debleach_polyfit_Callback(hObject, eventdata, handles)
% hObject    handle to debleach_polyfit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ord=str2num(get(handles.order,'String'));
st=get(handles.List_trace,'String');
val=get(handles.List_trace,'Value');
[str,remain]=strtok(st(val,:),'-');
[str,remain]=strtok(remain,'-');
%Border_Num =str2num(remain(2));
Border_Num=str2num(str{1,1});
%Tri = get(handles.Trial,'String');
%kk=str2double(Tri);
kk = get(handles.Trial_List,'Value');
Fr_Num = get(handles.Farme_Num,'String');
handles.current_s=handles.trace{2,val};
if ord==0
    plot(handles.axes2,handles.current_s);
   % figure,plot(handles.current_s)
else

[xx,p] = vscope_debleach(handles.current_s, ord);

plot(handles.axes2,xx);
set(handles.P,'String',p)
end

%    handles=handles.current_s;
%    save handles handles
%figure,plot(xx)
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Trial_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Trial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Save_P.
function Save_P_Callback(hObject, eventdata, handles)
% hObject    handle to Save_P (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName]=uiputfile(strcat('P\',datestr(now,1),num2str(hour(now)),num2str(minute(now)),'.mat'));
P1=get(handles.P,'String');
P=str2num(P1);
%borders=handles.borders;
%BWs=handles.BW;
save(strcat(PathName,FileName),'P'); 
guidata(hObject,handles)


% --- Executes on selection change in List_trace.
function List_trace_Callback(hObject, eventdata, handles)
% hObject    handle to List_trace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
st=get(hObject,'String');
val=get(hObject,'Value');
[str,remain]=strtok(st(val,:),'-');
[str,remain]=strtok(remain,'-');
%Border_Num =str2num(remain(2));
Border_Num=str2num(str{1,1});
%Tri = get(handles.Trial,'String');
%kk=str2double(Tri);
kk = get(handles.Trial_List,'Value');
Fr_Num = get(handles.Farme_Num,'String');
axes(handles.axes1)
imagesc(handles.A{kk}(:,:,str2num(Fr_Num))),colormap gray
hold on
%for i=1:size(Border_Num,2)
b=Border_Num(1);
boundary=handles.borders{1,b};
plot(boundary(:,2),boundary(:,1),'g','LineWidth',2)
%text(position(1,1),position(1,2),['\fontsize{12}\color{green}\leftarrow ' int2str(Border_Num)],...
%    'HorizontalAlignment','left')
text(boundary(1,2),boundary(1,1),['\fontsize{12}\color{green}' int2str(b) '\rightarrow '],...
    'HorizontalAlignment','right')
%end
hold off
plot(handles.axes2,handles.trace{2,val(1,1)}) 
%handles.current_s=s;
guidata(hObject,handles)
% Hints: contents = get(hObject,'String') returns List_trace contents as cell array
%        contents{get(hObject,'Value')} returns selected item from List_trace
%contents = get(hObject,'String');
%save contents
%adding new Member: con=[contents;'me'];

% --- Executes during object creation, after setting all properties.
function List_trace_CreateFcn(hObject, eventdata, handles)
% hObject    handle to List_trace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Plot_traces.
function Plot_traces_Callback(hObject, eventdata, handles)
% hObject    handle to Plot_traces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
List_borders =get(handles.border,'Value');%-1;
%Border_Num=str2double(Border_Nu);
%Tri = get(handles.Trial,'String');
%kk=str2double(Tri);
st=get(handles.Trial_List,'String');
kk = get(handles.Trial_List,'Value');
axes(handles.axes2)
ccd1=handles.A{kk};
%axes(handles.axes2)
%hold on
for i=1:size(List_borders,2)
    s=trace(ccd1,handles.BW{List_borders(i)});
    handles.trace_Num=handles.trace_Num+1;
    Name_s{1,1}=strcat(st{kk},'_Border','-',num2str(List_borders(i)));
    handles.trace{1,handles.trace_Num}=Name_s;% name trace inja zakhire mishavad
    handles.trace{2,handles.trace_Num}=s;% khode trace inja zakhire mishavad
    st2=get(handles.List_trace,'String');
    set(handles.List_trace,'String',[st2;Name_s]);
    %plot(s)
end
%hold off
%axes(handles.axes1)
guidata(hObject,handles)

% --- Executes on button press in delet_traces.
function delet_traces_Callback(hObject, eventdata, handles)
% hObject    handle to delet_traces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
trace_Nu=get(handles.List_trace,'Value');%-1;
trace_Numb=get(handles.List_trace,'String');% masalan 3dar7 char
kk=get(handles.Trial_List,'Value');
Fr_Num = get(handles.Farme_Num,'String');
%axes(handles.axes1)
%borders=handles.borders;
if trace_Nu(1,1)==1
    handles.trace=handles.trace(:,2:end);
    set(handles.List_trace,'String',trace_Numb(2:end,:),'Value',1);
else if trace_Nu(1,1)==size(trace_Numb,1)
        handles.trace=handles.trace(:,1:trace_Nu(1,1)-1);
        set(handles.List_trace,'String',trace_Numb(1:end-1,:),'Value',trace_Nu(1,1)-1);
    else 
   %     s=handles.borders;
   %     save s s
       handles.trace=[handles.trace(:,1:trace_Nu(1,1)-1),handles.trace(:,trace_Nu(1,1)+1:end)];
       set(handles.List_trace,'String',[trace_Numb(1:trace_Nu(1,1)-1,:);trace_Numb(trace_Nu(1,1)+1:end,:)],'Value',trace_Nu(1,1)-1);
    end
end
%Border_Numb=
%handles.border =set(handles.border,'String',);
%k=get(handles.border,'String');
%set(handles.border,'String',[k;Border_Num])
trace_Num=handles.trace_Num;
if trace_Num~=0
handles.trace_Num=handles.trace_Num-1;
end
%set(handles.List_trace,'String',[1:Border_Num]','Value',1)
guidata(hObject,handles)


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);

% --- Executes on button press in Save_All.
function Save_All_Callback(hObject, eventdata, handles)
% hObject    handle to Save_All (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName]=uiputfile(strcat('Save\',datestr(now,1),num2str(hour(now)),num2str(minute(now)),'.mat'));
borders=handles.borders;
BWs=handles.BW;
traces=handles.trace;
Trials(1,:)=handles.A;
NameTrials=get(handles.Trial_List,'String');
Trials(2,:)=NameTrials';
save(strcat(PathName,FileName),'borders','BWs','traces','Trials'); 
guidata(hObject,handles)

% --- Executes on button press in delet_Trial.
function delet_Trial_Callback(hObject, eventdata, handles)
% hObject    handle to delet_Trial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Trial_Nu=get(handles.Trial_List,'Value');%-1;
Trial_Numb=get(handles.Trial_List,'String');% masalan 3dar7 char
kk=get(handles.Trial_List,'Value');
Fr_Num = get(handles.Farme_Num,'String');
%axes(handles.axes1)
%borders=handles.borders;
if Trial_Nu(1,1)==1
    handles.A=handles.A(:,2:end);
    set(handles.Trial_List,'String',Trial_Numb(2:end,:),'Value',1);
else if Trial_Nu(1,1)==size(Trial_Numb,1)
        handles.A=handles.A(:,1:Trial_Nu(1,1)-1);
        set(handles.Trial_List,'String',Trial_Numb(1:end-1,:),'Value',Trial_Nu(1,1)-1);
    else 
   %     s=handles.borders;
   %     save s s
       handles.A=[handles.A(:,1:Trial_Nu(1,1)-1),handles.A(:,Trial_Nu(1,1)+1:end)];
       set(handles.Trial_List,'String',[Trial_Numb(1:Trial_Nu(1,1)-1,:);Trial_Numb(Trial_Nu(1,1)+1:end,:)],'Value',Trial_Nu(1,1)-1);
    end
end
%Border_Numb=
%handles.border =set(handles.border,'String',);
%k=get(handles.border,'String');
%set(handles.border,'String',[k;Border_Num])
handles.Trial_Num=handles.Trial_Num-1;
Trial_Num=handles.Trial_Num;
%set(handles.Trial_List,'String',[1:Border_Num]','Value',1)
guidata(hObject,handles)


% --- Executes on button press in plot_toghether.
function plot_toghether_Callback(hObject, eventdata, handles)
% hObject    handle to plot_toghether (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% to plot normalized signals toghether with the border related to them in
% colored way
color=[1,0,0;0,1,0;0,0,1;1,1,0;0,1,1;1,0,1;.5,0,0;0,.5,0;0,0,.5;.5,.5,0];
%[0.5,.5,1;0,0,1;0,1,1;0,1,0;1,0,0;1,0,1;1,1,0;1,1,.5;.5,.5,1;0,1,.5];
st=get(handles.List_trace,'String');
val=get(handles.List_trace,'Value');
[str,remain]=strtok(st(val,:),'-');
[str,remain]=strtok(remain,'-');
%Border_Num =str2num(remain(2));
%Border_Num=str2num(str{:});
Border_Num=str2num(cell2mat(str));
%Tri = get(handles.Trial,'String');
%kk=str2double(Tri);
kk = get(handles.Trial_List,'Value');
Fr_Num = get(handles.Farme_Num,'String');
axes(handles.axes1)
imagesc(handles.A{kk}(:,:,str2num(Fr_Num))),colormap gray

for i=1:size(Border_Num)
b=Border_Num(i);
boundary=handles.borders{1,b};
hold on
plot(boundary(:,2),boundary(:,1),'color',color(b,:),'LineWidth',2)
end
hold off

T=[];
axes(handles.axes2)

nf1=get(handles.nf,'string');% number of frames used for normalization
nf=str2num(nf1);
for i=1:size(val,2)
    b=val(i);
    traces=handles.trace{2,b};
    T=(traces-mean(traces(1:nf)))/std(traces(1:nf));
    plot(handles.axes2,T,'color',color(b,:))
    hold on
end
hold off
axes(handles.axes1)
%handles.current_s=s;
guidata(hObject,handles)



function nf_Callback(hObject, eventdata, handles)
% hObject    handle to nf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nf as text
%        str2double(get(hObject,'String')) returns contents of nf as a double


% --- Executes during object creation, after setting all properties.
function nf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to nf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nf as text
%        str2double(get(hObject,'String')) returns contents of nf as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in debleach.
function debleach_Callback(hObject, eventdata, handles)
% hObject    handle to debleach (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in deb_typ.
function deb_typ_Callback(hObject, eventdata, handles)
% hObject    handle to deb_typ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns deb_typ contents as cell array
%        contents{get(hObject,'Value')} returns selected item from deb_typ


% --- Executes during object creation, after setting all properties.
function deb_typ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to deb_typ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plot_togheter_debleached.
function plot_togheter_debleached_Callback(hObject, eventdata, handles)
% hObject    handle to plot_togheter_debleached (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
color=[1,0,0;0,1,0;0,0,1;1,1,0;0,1,1;1,0,1;.5,0,0;0,.5,0;0,0,.5;.5,.5,0];
%[0.5,.5,1;0,0,1;0,1,1;0,1,0;1,0,0;1,0,1;1,1,0;1,1,.5;.5,.5,1;0,1,.5];
st=get(handles.List_trace,'String');
val=get(handles.List_trace,'Value');
[str,remain]=strtok(st(val,:),'-');
[str,remain]=strtok(remain,'-');
%Border_Num =str2num(remain(2));
%Border_Num=str2num(str{:});
Border_Num=str2num(cell2mat(str));
%Tri = get(handles.Trial,'String');
%kk=str2double(Tri);
kk = get(handles.Trial_List,'Value');
Fr_Num = get(handles.Farme_Num,'String');
axes(handles.axes1)
imagesc(handles.A{kk}(:,:,str2num(Fr_Num))),colormap gray

for i=1:size(Border_Num)
b=Border_Num(i);
boundary=handles.borders{1,b};
hold on
plot(boundary(:,2),boundary(:,1),'color',color(b,:),'LineWidth',2)
end
hold off

T=[];
axes(handles.axes2)

nf1=get(handles.nf,'string');% number of frames used for normalization
nf=str2num(nf1);
for i=1:size(val,2)
    b=val(i);
    traces=handles.trace{2,b};
    T=(traces-mean(traces(1:nf)))/std(traces(1:nf));
    %*********************************
    y=T+20;% takhmine debleach bar ruye kole dadeha
    edata(:,1)=1:size(y,1);
%cfun = fit(edata(:,1),(y(:,1))-min(y)+1,'exp1');
    cfun = fit(edata(:,1),(y(:,1)),'exp1');
%cf2=polyfit(edata(:,1),y(:,1),1);
    y_deb = feval(cfun,edata(:,1))-20;
    T_deb=T-y_deb;
%*************************
    plot(handles.axes2,T_deb,'color',color(b,:))
    hold on
end
hold off
axes(handles.axes1)
%handles.current_s=s;
guidata(hObject,handles)



function Analog_Num_Callback(hObject, eventdata, handles)
% hObject    handle to Analog_Num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Analog_Num = get(hObject,'String');
Trial_Num=get(handles.Trial_List,'Value');
if ~isempty(handles.Analog{Trial_Num})
plot(handles.axes_Analog,handles.Analog{Trial_Num}(:,str2num(Analog_Num)))
end
% Hints: get(hObject,'String') returns contents of Analog_Num as text
%        str2double(get(hObject,'String')) returns contents of Analog_Num as a double


% --- Executes during object creation, after setting all properties.
function Analog_Num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Analog_Num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function varargout = JGFE_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function path_Callback(hObject, eventdata, handles)
% hObject    handle to path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of path as text
%        str2double(get(hObject,'String')) returns contents of path as a double


% --- Executes during object creation, after setting all properties.
function path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Analog_Channels_Callback(hObject, eventdata, handles)
% hObject    handle to Analog_Channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Analog_Channels as text
%        str2double(get(hObject,'String')) returns contents of Analog_Channels as a double


% --- Executes during object creation, after setting all properties.
function Analog_Channels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Analog_Channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
