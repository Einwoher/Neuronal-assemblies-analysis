 %the number of the trial
choice = menu('Choose a mode','Spike detection','Raster plot analysis');


currentFolder = pwd;
SaveRegion=fullfile(pwd,'region ');

listing1 = dir('*signals.mat');
load(char(listing1.name));
listing2 = dir('*regions.mat');
load(char(listing2.name));
   
[ii jj] = ndgrid(1:nx,1:ny);
regioncenters = A * [ii(:) jj(:)]; %the XY coordinates of all ROIs
% A=exist('region','var');

if choice==1 % for the mode Spike detection
    prompt = {'Enter a trial number'};
dlgtitle = 'Trial number';
definput = {'1'};
dims = [1 40];
opts.Interpreter = 'tex';
answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
TrNum = str2num(answer{1});
%     signals=NeuropilCorr(signals,localneuropil);
%     TrNum=1; %trial number
%     signals=NeuropilCorrLSB(signals,localneuropil,TrNum);%least square baseline correction of the signal
    signals=NeuropilCorr(signals,localneuropil); 
region.traces=transpose(signals(:,:,TrNum)); % matrix
region.fits=zeros(length(signals(1,:,TrNum)),length(signals(:,1,TrNum)));
region.name=cell(1,1);
region.contours=cell(1,length(regioncenters));
region.location=ones(1,length(regioncenters));
region.name = {'Name 1'};
region.contours=cell(1,length(regioncenters));
region.onsets=cell(1,length(regioncenters));  %cells
region.offsets=cell(1,length(regioncenters)); % cells



% dt=300/9543; %(time between frames)
% region.traces=transpose(ans.res.fit);

for i=1:length(regioncenters) %generation of circle contours arround the cell centers
    th = 0;
    for j=1:100
        th = th+pi/50;
    region.contours{i}(j,1)=5 * cos(th) + regioncenters(i,1);  %X of circle contour
    region.contours{i}(j,2)= 5 * sin(th) + regioncenters(i,2); %Y of circle contour
    end;
end;
end;
 if choice==2 %for the mode Raster plot analysis
     
     filename = uigetfile('*.mat');
     load(filename);
     filename = filename(1:end-4);
    
     
%     % randomise cell indexing to avoid clustering bias
%     indexvector=randperm(length(region.onsets));
%     region.onsets=region.onsets(indexvector); %ok
%     region.offsets=region.offsets(indexvector); %ok
%     region.fits=region.fits(indexvector,:);
%     region.contours=region.contours(indexvector); %ok
%     region.traces=region.traces(indexvector,:);
%     region.location=region.location(indexvector);
%     
 end;
 listing3 = dir('*avgimg.png');
importedData=uiimport(char(listing3.name));

   while ~exist('importedData','var')
      pause(1)
      fprintf('Waiting for importedData\n')
    end
region.image = cell2mat(struct2cell(importedData));
%% 
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
sz = get(0,'screensize');
tr = region.traces; %array of raw traces for all the cells in the field: Ncells x 999
nt = [];
for c = 1:size(tr,1)
%       nt(c,:) = tr(c,:); %in case of least square correction baseline
     nt(c,:) = dfoverf(tr(c,:))*100; %array of traces in df/f percentage Ncells x 999 
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

if choice==2
    TrNum=sscanf(filename,'region%d');
for c = 1:size(spk,1)
    spk(c,region.onsets{c}) = 1;
    dec(c,region.offsets{c}) = 1;
end
end;

%general window.
fig = figure('Name','SpikesAnalysis','NumberTitle','off','toolbar','figure','position',[1 0.15*sz(4) sz(3) 0.7*sz(4)],'doublebuffer','on'); %create a window to display images

%field of view window with contrast adjustment tools
uicontrol('Style','text','Units','normalized','String','Image','Position',[.02 .8 .11 0.03],'FontSize',12,'FontWeight','Bold',...
    'BackgroundColor',[.8 .8 .8]);
uicontrol('Style','text','units','normalized','string','Brightness','position',[.02 .75 .11 .02],'FontSize',9,'BackgroundColor',[.8 .8 .8]);
bbright = uicontrol('Style','slider','Units','normalized','Position',[.02 .72 .11 .02],'Min',0,'Max',1,'Sliderstep',[.01 .05],'Value',1/3, ...
    'Enable','off','Callback','HippoContrast');
uicontrol('Style','text','units','normalized','string','Contrast','position',[.02 .70 .11 .02],'FontSize',9,'BackgroundColor',[.8 .8 .8]);
bcontrast = uicontrol('Style','slider','Units','normalized','Position',[.02 .68 .11 .02],'Min',0,'Max',1,'Sliderstep',[.01 .05],'Value',1/3, ...
    'Enable','off','Callback','HippoContrast');

a=region.image;
imgax = subplot('position',[0.1 0.703571 0.203125 0.282143]); %Create axes in tiled positions
imagesc(a); %Scale data and display image object
hold on
set(gca,'xtick',[],'ytick',[]);
axis equal
axis tight
box on
colormap gray
set(bbright,'enable','on');
set(bcontrast,'enable','on');
[maxy maxx] = size(a);
% zoom on;
HippoContrast;

%buttons
plotonset = uicontrol('Style','pushbutton','Units','normalized','String','Plot onsets','Position',[.3 .7 .1 0.04],'FontSize',12,...
     'Callback','PlotOnsetsUnic');

Clustering = uicontrol('Style','pushbutton','Units','normalized','String','Clustering','Position',[.3 .8 .1 0.04],'FontSize',12,...
     'Callback','ClusteringAnton');
SpikeDetecting = uicontrol('Style','pushbutton','Units','normalized','String','All neurons spikes','Position',[.3 .9 .1 0.04],'FontSize',12,...
     'Callback','MLSpikeDetection');

 SpikeDetecting = uicontrol('Style','pushbutton','Units','normalized','String','One neuron spikes','Position',[.3 .85 .1 0.04],'FontSize',12,...
     'Callback','SingleCellSpikeDetection');
 
 plotonset = uicontrol('Style','pushbutton','Units','normalized','String','Propagation','Position',[0.7 .7 .1 0.04],'FontSize',12,...
     'Callback','Propagation');
 plotonset = uicontrol('Style','pushbutton','Units','normalized','String','Propagation','Position',[0.7 .8 .1 0.04],'FontSize',12,...
     'Callback','EventsPlot');
 
 uicontrol('Style','text','Units','normalized','String','Event #','Position',[.79 0.75 .05 0.04],'FontSize',12,'FontWeight','Bold',...
    'HorizontalAlignment','right','BackgroundColor',[.8 .8 .8]);
uicontrol('Style','text','Units','normalized','String','From','Position',[.79 0.85 .05 0.04],'FontSize',12,'FontWeight','Bold',...
    'HorizontalAlignment','right','BackgroundColor',[.8 .8 .8]);
uicontrol('Style','text','Units','normalized','String','To','Position',[.84 0.85 .05 0.04],'FontSize',12,'FontWeight','Bold',...
    'HorizontalAlignment','right','BackgroundColor',[.8 .8 .8]);
 EventsToPlot = uicontrol('Style','edit','Units','normalized','String','all','Position',[0.8 .7 .05 0.04],'FontSize',12,'FontWeight','Bold',...
      'BackgroundColor',[1 1 1],'HorizontalAlignment','left');
  EventStart = uicontrol('Style','edit','Units','normalized','String','0','Position',[0.8 .8 .05 0.04],'FontSize',12,'FontWeight','Bold',...
      'BackgroundColor',[1 1 1],'HorizontalAlignment','left');
  EventEnd = uicontrol('Style','edit','Units','normalized','String','0','Position',[0.86 .8 .05 0.04],'FontSize',12,'FontWeight','Bold',...
      'BackgroundColor',[1 1 1],'HorizontalAlignment','left');
  
 % visualisation of traces
 xlimits = [0 size(nt,2)+1];
 ylimits = [mingl maxgl];



positionVector1 = [0.03, 0.50, 0.94, 0.17];    % position of first subplot
box on
set(gca,'buttondownfcn','hevZoom')
% 
trax=subplot('Position',positionVector1);

%plot(x,y1) aus.res.spikes


positionVector2 = [0.03, 0.28, 0.94, 0.17];    % position of second subplot
subplot('Position',positionVector2)
 image(transpose(signals(:,:,TrNum)),'CDataMapping','scaled')
%plot(x,y1)

positionVector3 = [0.03, 0.06, 0.94, 0.17];    % position of second subplot
subplot('Position',positionVector3)
% plot(ans.res.fit);
% Plot(ans.res.spikest/dt)


% trax = subplot('position',[0.03 0.2 0.94 0.480357]);
% % box on
% % set(gca,'buttondownfcn','hevZoom')
% 
% fitplot=subplot('position',[0.03 0.1 0.94 0.480357]);
% % box on
% % set(gca,'buttondownfcn','hevZoom')

def=0;
allautomatic=0;

onsetmat(1:size(nt,1),1:length(nt(1,:)))=0;
offsetmat(1:size(nt,1),1:length(nt(1,:)))=0;
uicontrol('Style','text','Units','normalized','String','Cell #','Position',[.05 0.95 .05 0.04],'FontSize',12,'FontWeight','Bold',...
    'HorizontalAlignment','right','BackgroundColor',[.8 .8 .8]);
txcellnum = uicontrol('Style','edit','Units','normalized','String','1','Position',[.05 .9 .05 0.04],'FontSize',12,'FontWeight','Bold',...
       'BackgroundColor',[1 1 1],'HorizontalAlignment','left','Callback','hevPlotTraceHpEvWithNetActUNIC');
 
bgoto = uicontrol('Style','pushbutton','Units','normalized','String','Go','Position',[.05 .85 .05 0.04],'FontSize',12,...
    'Callback','hevPlotTraceHpEvWithNetActUNIC'); 


xlimits = [0 size(nt,2)+1];
ylimits = [mingl maxgl];
hevPlotTraceHpEvWithNetActUNIC %automatically plot the traces for the first cell

subplot(imgax) %first cell countour plot
cellSel=num;

for i= cellSel
    hold on;
    plot(region.contours{i}([1:end 1],1),region.contours{i}([1:end 1],2),'r-')
end
%% ML spike parameters
uicontrol('Style','text','Units','normalized','String','tau','Position',[.44 0.9 .05 0.04],'FontSize',12,'FontWeight','Bold',...
    'HorizontalAlignment','right','BackgroundColor',[.8 .8 .8]);
tau = uicontrol('Style','edit','Units','normalized','String','1.87','Position',[0.5 .9 .05 0.04],'FontSize',12,'FontWeight','Bold',...
      'BackgroundColor',[1 1 1],'HorizontalAlignment','left');
  uicontrol('Style','text','Units','normalized','String','drift','Position',[.44 0.85 .05 0.04],'FontSize',12,'FontWeight','Bold',...
    'HorizontalAlignment','right','BackgroundColor',[.8 .8 .8]);
drift0 = uicontrol('Style','edit','Units','normalized','String','10','Position',[0.5 .85 .05 0.04],'FontSize',12,'FontWeight','Bold',...
      'BackgroundColor',[1 1 1],'HorizontalAlignment','left');
  
% uicontrol('Style','text','Units','normalized','String','Saturation','Position',[.44 0.80 .05 0.04],'FontSize',12,'FontWeight','Bold',...
%     'HorizontalAlignment','right','BackgroundColor',[.8 .8 .8]);
% saturation = uicontrol('Style','edit','Units','normalized','String','0.1','Position',[0.5 .80 .05 0.04],'FontSize',12,'FontWeight','Bold',...
%       'BackgroundColor',[1 1 1],'HorizontalAlignment','left');  
%   
%   
%   uicontrol('Style','text','Units','normalized','String','A','Position',[.44 0.75 .05 0.04],'FontSize',12,'FontWeight','Bold',...
%     'HorizontalAlignment','right','BackgroundColor',[.8 .8 .8]);
% a = uicontrol('Style','edit','Units','normalized','String','0.1','Position',[0.5 .75 .05 0.04],'FontSize',12,'FontWeight','Bold',...
%       'BackgroundColor',[1 1 1],'HorizontalAlignment','left');  

  