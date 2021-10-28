% clear;
% [filename, pathname] = uigetfile({'*.mat'}, 'Choose file to open');
% if ~isstr(filename)
%     return
% end
% fnm = [pathname filename];
% load (fnm)
[ii jj] = ndgrid(1:nx,1:ny);
regioncenters = A * [ii(:) jj(:)]; %the XY coordinates of all ROIs
region.traces=transpose(signals(:,:,1)); % matrix
 xStart = 1;
   dx = 1;
   N = length(signals(1,:,:));
region.userdata.active=xStart + (0:N-1)*dx;
region.userdata.neurons=xStart + (0:N-1)*dx;
% to run thomas spike identification algorithm size(signals) times and to
% create the matrix of onsets for all cells taking ans.res.spikes generated
% by this algorithm. 
%region.traces=transpose(ans.res.fit);

region.contours=cell(1,length(regioncenters));

for i=1:length(regioncenters) %generation of circle contours arround the cell centers
    th = 0;
    for j=1:100
        th = th+pi/50;
    region.contours{i}(j,1)=5 * cos(th) + regioncenters(i,1);  %X of circle contour
    region.contours{i}(j,2)= 5 * sin(th) + regioncenters(i,2); %Y of circle contour
    end;
end;
% region.image = x1710100x2Dreg0x2Davgimg;
% region.image = x171206_0x2Dreg0x2Davgimg
% region.image = x16012018_0x2Dreg0x2Davgimg
% region.image = x24012018_0x2Dreg0x2Davgimg;
%region.image = x13022018_mouse0x2Dreg0x2Davgimg;
% region.image = x070218_mouse0x2Dreg0x2Davgimg
% region.image =x13022018_mouse0x2Dreg0x2Davgimg
% region.image =x220218_m0x2Dreg0x2Davgimg
region.image =x180307_m0x2Dreg0x2Davgimg;
%reregion.image =x180307_m0x2Dreg0x2Davgimg;gion.location=regioncenters;
%viscircles(regioncenters,5)
%region.image = x171010_Reg_Avgimg;
% region.image = x1710100x2Dreg0x2Davgimg; % av image
% region.image = x171201A1_Reg_Avgimg;
% region.image = x171206_0x2Dreg0x2Davgimg;
% region.onsets=cell(1,length(regioncenters));  %cells
% region.offsets=cell(1,length(regioncenters)); % cells
plateauindex=0;
schmutz=0;
indsave=0;
indsave1=0;
defnext=0;
ind0=0;
schm=0;
numOld=1;
shi=0;
saveind=0;
yind=0;
% opengl neverselect; %disables autoselection of OpenGL.
sz = get(0,'screensize');
fig = figure('Name','HippoEvents','NumberTitle','off','toolbar','figure','position',[1 0.15*sz(4) sz(3) 0.7*sz(4)],'doublebuffer','on'); %create a window to display images

% [filename, pathname] = uigetfile({'*.mat'}, 'Choose data file to open')
% 
% if ~isstr(filename)
%     delete(gcf)
%     clear
% end
% fnm = [pathname filename];
% load(fnm)

%gray scale image adjustment interface
uicontrol('Style','text','Units','normalized','String','Image','Position',[.87 .955 .11 0.03],'FontSize',12,'FontWeight','Bold',...
    'BackgroundColor',[.8 .8 .8]);
uicontrol('Style','text','units','normalized','string','Brightness','position',[.87 .88 .11 .02],'FontSize',9,'BackgroundColor',[.8 .8 .8]);
bbright = uicontrol('Style','slider','Units','normalized','Position',[.87 .86 .11 .02],'Min',0,'Max',1,'Sliderstep',[.01 .05],'Value',1/3, ...
    'Enable','off','Callback','HippoContrast');
uicontrol('Style','text','units','normalized','string','Contrast','position',[.87 .83 .11 .02],'FontSize',9,'BackgroundColor',[.8 .8 .8]);
bcontrast = uicontrol('Style','slider','Units','normalized','Position',[.87 .81 .11 .02],'Min',0,'Max',1,'Sliderstep',[.01 .05],'Value',1/3, ...
    'Enable','off','Callback','HippoContrast');
nextident = uicontrol('Style','pushbutton','Units','normalized','String','Identify next','Position',[.63 .92 .11 0.04],'FontSize',12,...
    'Callback', 'def = 1; hevButtonDownNetAct');
schmutzconfirm = uicontrol('Style','pushbutton','Units','normalized','String','Confirm schmutz','Position',[.63 .87 .11 0.04],'FontSize',12,...
    'Callback', 'schm = 1; hevPlotTraceHpEvWithNetActAnton');
 manual = uicontrol('Style','pushbutton','Units','normalized','String','Manual analysis','Position',[.63 .77 .11 0.04],'FontSize',12,...
     'Callback', 'manual_analysis; button');
 set(manual,'enable','off');

a=region.image;
%a=scatter(regioncenters(:,1),regioncenters(:,2));%region.image;

imgax = subplot('position',[0.4234375 0.703571 0.203125 0.282143]); %Create axes in tiled positions
imagesc(a); %Scale data and display image object
hold on

set(gca,'xtick',[],'ytick',[]);
axis equal
axis tight
box on

colormap gray

% set(bzoom,'enable','on');
set(bbright,'enable','on');
set(bcontrast,'enable','on');


[maxy maxx] = size(a);

% zoom on;
HippoContrast;

%traces
%schmutz=region.userdata.schmutzr;

tr = region.traces; %array of raw traces for all the cells in the field: Ncells x 999
nt = [];
for c = 1:size(tr,1)
    nt(c,:) = -1*dfoverf(tr(c,:))*100; %array of traces in df/f percentage Ncells x 999 %multiplied by -1 to be in negative
    minmax(c,1) = min(nt(c,:)); %array of the max and min df/f for all the cells in the field
    minmax(c,2) = max(nt(c,:));
end
mingl = min(minmax(:,1)); %the lowest df/f in the field
if max(minmax(:,2))<=200
maxgl = max(minmax(:,2)); %the highest df/f in the field
end;

if max(minmax(:,2))>200
maxgl = 200; %the highest df/f in the field
end;

spk = zeros(size(nt));
dec = zeros(size(nt));
for c = 1:size(spk,1)
    spk(c,region.onsets{c}) = 1;
    dec(c,region.offsets{c}) = 1;
end

for i=1:length(nt(:,1))
if isempty(region.onsets{i})==0
    plateauindex=1;
    break
end;
end

% if  plateauindex==1
% schmutz=region.userdata.schmutzr;
% shi=length(region.userdata.schmutzr);
% plateauindex=0;
% end;


xlimits = [0 size(nt,2)+1];
ylimits = [mingl maxgl];

trax = subplot('position',[0.03 0.2 0.94 0.480357]);
box on
set(gca,'buttondownfcn','hevZoom')
%cd('C:\Program Files\MATLAB\R2006b\work\Prog_Anton');
cd('C:\Users\Anton\Documents\MATLAB\work\Prog_Anton');

%cd('\\Netdata\Rotarien\CARRON Romain\work_ANTON\R2006b\work\Prog_Anton');
%baseline interface:
def=0;
allautomatic=0;

onsetmat(1:size(nt,1),1:length(nt(1,:)))=0;
offsetmat(1:size(nt,1),1:length(nt(1,:)))=0;

uicontrol('Style','text','Units','normalized','String','Baseline adjustment ','Position',[.05 .955 0.11 0.03],'FontSize',12,'FontWeight','Bold',...
    'HorizontalAlignment','right','BackgroundColor',[.8 .8 .8]);
 uicontrol('Style','text','units','normalized','string','P','position',[.05 .835 .11 .02],'FontSize',9,'BackgroundColor',[.8 .8 .8]);
 padj = uicontrol('Style','slider','Units','normalized','Position',[.05 .81 .11 .02],'Min',0.001,'Max',0.9999,'Sliderstep',[.001 .01],'Value',0.6, ...
     'Enable','off','Callback','ind0=1; hevPlotTraceHpEvWithNetActAnton');
  uicontrol('Style','text','units','normalized','string','Threshold','position',[.05 .935 .11 .02],'FontSize',9,'BackgroundColor',[.8 .8 .8]);
  threshold = uicontrol('Style','slider','Units','normalized','Position',[.05 .91 .11 .02],'Min',0.00001,'Max',1,'Sliderstep',[.001 .01],'Value',0.5, ...
     'Enable','off','Callback','ind0=1; hevPlotTraceHpEvWithNetActAnton');

 uicontrol('Style','text','units','normalized','string','Lambda','position',[.05 .785 .11 .02],'FontSize',9,'BackgroundColor',[.8 .8 .8]);
lambdaadj = uicontrol('Style','slider','Units','normalized','Position',[.05 .76 .11 .02],'Min',0.0000000001,'Max',1,'Sliderstep',[0.00001 0.0001],'Value',0.5, ...
    'Enable','off','Callback','ind0=1; hevPlotTraceHpEvWithNetActAnton');
uicontrol('Style','text','units','normalized','string','Event duration','position',[.05 .885 .11 .02],'FontSize',9,'BackgroundColor',[.8 .8 .8]);
eventduradj = uicontrol('Style','slider','Units','normalized','Position',[.05 .86 .11 .02],'Min',0.01,'Max',0.5,'Sliderstep',[0.01 0.1],'Value',0.3, ...
    'Enable','off','Callback','ind0=1; hevPlotTraceHpEvWithNetActAnton');
set(padj,'enable','on');
set(lambdaadj,'enable','on');
set(threshold,'enable','on');
set(eventduradj,'enable','on');

default = uicontrol('Style','pushbutton','Units','normalized','String','Default set','Position',[.05 .70 .11 0.04],'FontSize',12,...
    'Callback', 'def = 1; hevPlotTraceHpEvWithNetActAnton');

startlimit = uicontrol('Style','slider','Units','normalized','Position',[.021 .15 0.95 .02],'Min',0,'Max',1,'Sliderstep',[.001 .01],'Value',0.002, ...
     'Enable','off','Callback','ind0=1; hevPlotTraceHpEvWithNetActAnton');
% automatic=uicontrol('Style','pushbutton','Units','normalized','String','Automatic event search','Position',[.25 .05 .14 0.04],'FontSize',12,...
%     'Callback', 'AutoEventIdentification'); 

plotonset = uicontrol('Style','pushbutton','Units','normalized','String','Plot onsets','Position',[.54 .05 .1 0.04],'FontSize',12,...
     'Callback','PlotOnset');
plotonset = uicontrol('Style','pushbutton','Units','normalized','String','Clustering','Position',[.44 .05 .1 0.04],'FontSize',12,...
     'Callback','ClusteringAnton');
%=================================
%===y axis adjustment
yaxis = uicontrol('Style','slider','Units','normalized','Position',[.005 .2 .01 .48],'Min',0.001,'Max',1,'Sliderstep',[.01 .1],'Value',1, ...
     'Enable','off','Callback','yind=1; hevPlotTraceHpEvWithNetActAnton');
 set(yaxis,'enable','on');
 set(startlimit,'enable','on');
%===================== 
 
uicontrol('Style','text','Units','normalized','String','Cell #','Position',[.05 .05 .05 0.04],'FontSize',12,'FontWeight','Bold',...
    'HorizontalAlignment','right','BackgroundColor',[.8 .8 .8]);

  txcellnum = uicontrol('Style','edit','Units','normalized','String','1','Position',[.11 .05 .05 0.04],'FontSize',12,'FontWeight','Bold',...
      'BackgroundColor',[1 1 1],'HorizontalAlignment','left','Callback','hevPlotTraceHpEvWithNetActAnton');

bgoto = uicontrol('Style','pushbutton','Units','normalized','String','Go','Position',[.17 .05 .05 0.04],'FontSize',12,...
    'Callback','hevPlotTraceHpEvWithNetActAnton'); 

progtx = uicontrol('Style','text','Units','normalized','String','','Position',[.70 .05 .25 0.04],'FontSize',12,'FontWeight','Bold',...
    'HorizontalAlignment','left','BackgroundColor',[.8 .8 .8]);

currdir = pwd;
cd('C:\Users\Anton\Documents\MATLAB\work\Hippo\SignalDetectors');

%cd('C:\Program Files\MATLAB\R2006b\work\Hippo\SignalDetectors');
mt = dir('*.m');
cd(currdir);

st = cell(1,length(mt));
for c = 1:length(mt)
    st{c} = mt(c).name(1:end-2);
    if strcmp(upper(st{c}(1:min([8 length(st{c})]))),'Event_identificator')
        st{c} = st{c}(9:end);
        
        
    end
end

dummy(2) = uicontrol('Style','text','Units','normalized','String','Trace reader','Position',[.80 0.08 .11 0.04],'FontSize',9,...
    'HorizontalAlignment','left','BackgroundColor',[.8 .8 .8]);
dpreaders = uicontrol('Style','popupmenu','Units','normalized','String',st,'Position',[.80 .06 .18 0.025],'FontSize',9,...
    'BackgroundColor',[1 1 1]);

% bdetect1 = uicontrol('Style','pushbutton','Units','normalized','String','Detect current usual','Position',[.25 .05 .14 0.04],'FontSize',12,...
%     'Callback','spk(num,:) = 0; dec(num,:) = 0; [s d] = hippodettrial(tr(num,:)); spk(num,s) = 1; dec(num,d) = 1; hevPlotTraceHpEvWithNetActAnton');
% bdetect2 = uicontrol('Style','pushbutton','Units','normalized','String','Detect current Filt','Position',[.40 .1 .14 0.04],'FontSize',12,...
%     'Callback','spk(num,:) = 0; dec(num,:) = 0; [s d] = HippoEvent_DetSingTrHP(region,num,''no''); spk(num,s) = 1; dec(num,d) = 1; hevPlotTraceHpEvWithNetActAnton;');
% bdetect3 = uicontrol('Style','pushbutton','Units','normalized','String','Detect current Fast&Slow','Position',[.40 .005 .14 0.04],'FontSize',12,...
%     'Callback','FastAndSlowDetectNew');
% bdetect4 = uicontrol('Style','pushbutton','Units','normalized','String','better noise estimation','Position',[.61 .1 .14 0.04],'FontSize',12,...
%     'Callback','spk(num,:) = 0; dec(num,:) = 0; [s d energ] = HippoEvent_NoiseEstFastBest(region,num,''no'',region.onsets{num},region.offsets{num}); spk(num,s) = 1; dec(num,d) = 1; hevPlotTraceHpEvWithNetActAnton');
% 
% %fast current detection
% bdetect5 = uicontrol('Style','pushbutton','Units','normalized','String','Detect current Fast','Position',[.61 .005 .14 0.04],'FontSize',12,...
%     'Callback','spk(num,:) = 0; dec(num,:) = 0; [s d] = HippoEvent_FastBest(region,num,''no'',region.onsets{num},region.offsets{num}); spk(num,s) = 1; dec(num,d) = 1; hevPlotTraceHpEvWithNetActAnton');

% handle=uicontrol('Style','pushbutton','Units','normalized','String','uniform df/f scale','Position',[.61 .05 .1 0.04],'FontSize',12,...
%     'Callback','hevPlotTraceHpEvWithNetActAnton'); 




%====================

% bdetect = uicontrol('Style','pushbutton','Units','normalized','String','Detect all','Position',[.80 .02 .07 0.04],'FontSize',12,...
%     'Callback','hippoDispTrialNetActData');
bdeleteall = uicontrol('Style','pushbutton','Units','normalized','String','Delete events','Position',[.34 .05 .1 0.04],'FontSize',12,...
    'Callback','spk(num,:) = 0; dec(num,:) = 0; hevPlotTraceHpEvWithNetActAnton');
bsave = uicontrol('Style','pushbutton','Units','normalized','String','Save all','Position',[.63 .82 .11 0.04],'FontSize',12,...
    'Callback','saveind=1; hevSaveWithNetActData');
corrpair = uicontrol('Style','pushbutton','Units','normalized','String','Correlation','Position',[.63 .72 .11 0.04],'FontSize',12,...
    'Callback','Correlation');
 set(corrpair,'enable','off');
% if saveind==1
% set(manual,'enable','on')
% saveind=0;
% end;
% if defnext==1
% set(textcellnum,'string', num+1);
% end;

hevPlotTraceHpEvWithNetActAnton %automatically plot the traces for the first cell

subplot(imgax) %first cell countour plot
cellSel=num;
for i=cellSel
    hold on;
    plot(region.contours{i}([1:end 1],1),region.contours{i}([1:end 1],2),'r-')
end


