
% T = tps_trial('Z:\RAW DATA\Anton\180509\elphy 2P syncronisation.mesc');
T = tps_trial('G:\BackUp ANTON\Anton\raw data\180509\elphy 2P syncronisation.mesc');

numres=100; sig=1; p_val=0.05; %numres should be 200, initially was 1000
 clear('s')
 clear('episodes_res')
  dur=10; %jitter window size 20
  Window1=1; %sliding window size 10
  EventDuration = 1; % the minimal duration of the event to be taken into account 15
all_s = cell(1,2);
all_s{1} = zeros(size(region.traces)); % all cells
for i = 1:size(region.traces,1)
    all_s{1}(i,region.onsets{i})=1;
end
  s = all_s{1};
  
  episodes_1=zeros(1, length(nt));
    s_thick = s; %gauss_events(s,sig); %generes gaussian probability at time moments corresponding to the experimental data
    testprob=s_thick;
  nRow=size(s_thick,1);
  nCol= size(s_thick,1);
  nCol1=size(s_thick',1);
    trans=s_thick';
  corrs = zeros(size(s_thick,1));
 

 
 
%% generation of artificial dataset
ArtSpTr=cell(1,size(s,1));
 episodes = s;
    corrs_cont = reshape(corrs,1,prod(size(corrs))); %array of number of cells x nunber of cells of the elements of corr matrix  = experimental correlation
   
    episodes_cont = sum(episodes);
    corrs_res = zeros(fix(p_val*numres)+1,prod(size(corrs))); %matrix of 51 x prod(size(corrs) more than 5 percent of 1000 trials
    episodes_res=zeros(fix(p_val*numres)+1, length(s));
%% 
for t = 1:numres
     AllValDif=0;
     s = all_s{1};
     
                       while AllValDif==0   
    
        for c = 1:size(s,1); %size(s,1) number of cells
             ff=find(s(c,:)==1);
                      ff0=ff;   
              if length(ff)>1
              intervals=diff(ff);
                intervals(length(intervals)+1)=ff(1);
              
%              intervals(length(intervals)+1)=(9543-ff(length(ff)));
%              ix = randperm(length(intervals)); % randomise the indices
            intervals =intervals(randperm(length(intervals)));
%                      ff0=round((ff-dur)+rand(1,1)*(dur+dur)); 
%                      if (min(ff0)<=0)
%                         ff0(numel(min(ff)))=round((ff(1)+1)+rand(1,1)*(dur+dur));
%                      end;
%                      if max(ff0)>=(size(region.traces,2))
%                          ff0(numel(max(ff)))=round((ff(length(ff))-dur-dur-1)+rand(1,1)*(dur+dur));
%                     end;
%               intshuffled = intervals(ix);
              s(c,:)=0;
              ArtSpTr{c}(t,1:size(s,2))=0;
%              intshuffled = cumsum(intshuffled); 
               intshuffled = cumsum(intervals); 
             s(c,intshuffled)=1;
              ArtSpTr{c}(t,intshuffled)=1;
%                        intshuffled=ff0;
%                      if length(unique(intshuffled))==length(ff)&&(max(ff0))<size(region.traces,2)&&(min(ff0)>0) %check if the generated values doesn't repeat
                          AllValDif=1; 
%                          s(c,ff0)=1;
%                      end        
              end;
              if length(ff)<=1
              s(c,:) = s(c,randperm(size(s,2))); %size(s,2) is 999 - number of temporal frames, it returns the array of onsets put in random order 
               ArtSpTr{c}(t,:)=s(c,randperm(size(s,2)));
                        AllValDif=1; 
              end;
        end
     
                        end %of while       
        
        
        s_thick = s; %gauss_events(s,sig); %generate probabilities for the random distribution of onsets for during 1000 trials 
      
      corrs = zeros(size(s_thick,1));     %matrix of correlation for the artificial data set   

  
        corrs_res(1,:) = reshape(corrs,1,numel(corrs));
        episodes_0=sum(s);
       
        for i=1:length(episodes_0)
          episodes_1(i)=0;  
          AvInd=0;
             if ((i-Window1)>0)&&((i+Window1)<length(episodes_0))
              for ii=(i-Window1):(i+Window1)
%              episodes_1(i)=episodes_0(i-1)+episodes_0(i)+episodes_0(i+1);
               episodes_1(i)=episodes_1(i)+episodes_0(ii);
              AvInd=AvInd+1;
              end;
             end;  
        
           if ((i-Window1)<=0)
             for ii=1:(i+Window1)          
              episodes_1(i)=episodes_1(i)+episodes_0(ii);
             AvInd=AvInd+1;
            end;   
           end;
          
           if ((i+Window1)>=length(episodes_0))
              for ii=(i-Window1):length(episodes_0)
              episodes_1(i)=episodes_1(i)+episodes_0(ii);
              AvInd=AvInd+1;
             end; 
           end;
               
          
            episodes_1(i)=(episodes_1(i)/AvInd);
           
        end;  
        episodes_res(t,:) =  episodes_1;
%          episodes_res = sort(episodes_res,1);
%         corrs_res = sort(corrs_res); %array of probabilistic coincidences for 1000 randomly generates cases
        if mod(t,1) == 0 %every rial a point is printed out
            fprintf('.');
             size(find(s==1))
        end
end
    fprintf('\n');
   
 episodes_res = sort(episodes_res,1);
%%     
% filtering of real spike trains     
      for i=1:length(episodes_cont)
          episodes_1(i)=0;  
          AvInd=0;
             if ((i-Window1)>0)&&((i+Window1)<length(episodes_cont))
              for ii=(i-Window1):(i+Window1)
%              episodes_1(i)=episodes_0(i-1)+episodes_0(i)+episodes_0(i+1);
               episodes_1(i)=episodes_1(i)+episodes_cont(ii);
              AvInd=AvInd+1;
              end;
             end;  
        
           if ((i-Window1)<=0)
             for ii=1:(i+Window1)          
              episodes_1(i)=episodes_1(i)+episodes_cont(ii);
             AvInd=AvInd+1;
            end;   
           end;
          
           if ((i+Window1)>=length(episodes_cont))
              for ii=(i-Window1):length(episodes_cont)
              episodes_1(i)=episodes_1(i)+episodes_cont(ii);
              AvInd=AvInd+1;
             end; 
           end;
               
          
            Episodes_RealFiltered(i)=(episodes_1(i)/AvInd);
        end;  
    
    
    
  %episodes_res(96,:) threshold line
%   plotonset = figure('position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
%   raster = subplot('position',[0.03 0.2 0.94 0.3]);
%   subplot(raster)
% plot(Episodes_RealFiltered, 'gr')  
% hold on
% plot(episodes_res(96,:),'r')
% hold off
%% baseline


 lambda=100000000;

    p=0.01;
%     thhold=5;
%     eventdur=3;
%     set(padj,'value',0.6);
%     set(lambdaadj,'value',0.5);
%      set(threshold,'value',0.5);
%      set(eventduradj,'value',0.3);
%        set(startlimit,'value',0.002);

basel=baseline_(Episodes_RealFiltered',lambda,p);
% basel=basel-min(basel);
%%

% % % 
%  episodes_res(99,:)= basel+mean(episodes_res(99,:));   
%  EventThreshold= episodes_res(99,:);
  SpontThreshold=mean(episodes_res(99,:));
%  save ('ThresholdISO','SpontThreshold');
% %% for detection of evoked responses that trigger population events using the spontaneous threshold
episodes_res(99,:)= basel+SpontThreshold;   
EventThreshold= episodes_res(99,:);
%% detecting the events that are above the  dynamic sliding window  threshold 
clear('allcells');
clear('EventsIdxEd');
clear('EventsIdx');
clear('AboveThEvents');
%% new method of peaks identification 1 Savitzky-Golay filtering 2 local maxima 3 local minima 4 filtering over threshold
clear('EvTime')
clear('EvStart')
clear('EvEnd')
clear('EvDur')


y = sgolayfilt(Episodes_RealFiltered,3,7);
[pksmax,locsmax] = findpeaks(y); % local maxima
[pksmin,locsmin] = findpeaks(-y); %local minima
if locsmax(1)>locsmin(1)
    stpoint=1;
    endpoint=length(pksmax)-1;
    a=0;
    b=1;
    
end;
if locsmax(1)<locsmin(1)
    stpoint=2;
    endpoint=length(pksmax);
    a=-1;
    b=0;
end;

if locsmax(length(locsmax))>locsmin(length(locsmin))
   endpoint=endpoint-1; 
end;
% EventThreshold
j=1;
for i= stpoint:endpoint
    if pksmax(i)>EventThreshold(locsmax(i))
     EvTime(j)= locsmax(i);
     EvStart(j)=locsmin(i+a);
     EvEnd(j)=locsmin(i+b);
     EvDur(j)= EvEnd(j)-EvStart(j);
     j=j+1;
    end;
end;


for i=1:length(EvTime)
    jj=1;
    for j=EvStart(i)+1:EvEnd(i)
        AboveThEvents{i}(jj)=j;
        jj=jj+1;
     
       
    end;
end;

%% 
%      ThEvents=mat2cell( ThPoints, 1, diff( [0, find(diff(ThPoints) ~= 1), length(ThPoints)] )); %separate the array of points of the threshold below the signal line into cells of consequitive values
%      AboveThEvents = mat2cell( AboveThPoints, 1, diff( [0, find(diff(AboveThPoints) ~= 1), length(AboveThPoints)] )) ; %separate the array of points above the Threshold into cells of consequitive values
%      AboveThEvents = mat2cell( AboveThPoints, 1, diff( [0, find(diff(AboveThPoints) ~= 1), length(AboveThPoints)] )) ;
     EventsIdx=find(cellfun('length',AboveThEvents)>EventDuration); %chose only events that are longer than 1 point
      
%    m=1; 
%    for i=1:length(EventsIdx)  
%        j=EventsIdx(i);
%    if  mean(Episodes_RealFiltered(AboveThEvents{j}))-mean(EventThreshold(ThEvents{j}))>0.05*mean(Episodes_RealFiltered(AboveThEvents{j})) %the second condition fot the event, difference with the baseline should be higher than 5%
%        EventsIdxEd(m)=j;
%        m=m+1;
%    end
%    end
%%    
   clear('CommonFraction');
    clear('ComFrSize');
  clear('allcells');
 EventsIdxEd=length(EvTime);   
for i=1:length(EvTime);
%     length(EventsIdxEd)
%     j=EventsIdxEd(i);
      allcells{i}=[];
%      CommonFraction{i}=[];
    m=1;
%     for k=min(AboveThEvents{j}):max(AboveThEvents{j})
   for k=EvStart(i):EvEnd(i) % the adjacent events might overlap because the end of the firts is the beginning of the next one
        for n=1:length(all_s{1}(:,k)) %scan the columns of all_s matrix, each of which is time moment over all the cells
            if all_s{1}(n,k)==1
             allcells{i}(m)=n;
             m=m+1;
            end
        end
        
    end
%     EvStart(i)=min(AboveThEvents{j}); %% time when the event starts
%     EvEnd(i)=max(AboveThEvents{j}); %%time when the event ends
%     EvDur(i)=EvEnd(i)-EvStart(i)+1; %%duration of each event in frames 
   
    allcells{i}=unique(allcells{i}); 
end; 
%artificial spike trains participation
%%
CellCount=ParticipationThresholdSp(numres,ArtSpTr,EvTime,EvEnd,EvStart);
% CellCountSp=ParticipationThresholdSp(numres,ArtSpTr(ParticipativeCells),EvTime,EvEnd,EvStart);




%% Range of neurons by their participation rate in population events
A=cell2mat(allcells); %vector array from allcells containing the neurons participating in each population event
TNeur = tabulate(A); %attribute to each neuron its participation rate
[~,idx] = sort(TNeur(:,2),'descend'); % sort the neurons according to their participation rate, from min to max
sortedmat = TNeur(idx,:); % sort the neurons according to their participation rate, from min to max
sortedmat(1); %the most participative neuron

% filtering of cells according to their participation in population events
j=1;
for i=1:length(sortedmat);%length(CellCount)
  
    if sortedmat(i,2,1)>CellCount(sortedmat(i,1,1))
%    NewSortedMat(j)=sortedmat(i,1,1);
   
     NewSortedMat(j,1)=sortedmat(i,1,1);
     NewSortedMat(j,2)=sortedmat(i,2,1);
    j=j+1;
    end
end;
WrongCells=setdiff(sortedmat(:,1),NewSortedMat(:,1)); %all the cells that do not pass the filter
% WrongCells=[]; 
WrongCellsArray=cell(1,length(allcells));
for i=1:length(allcells)
  WrongCellsArray{i}=WrongCells;
end;
 
WrongCellsArrayPrctl=cell(1,5);
for i=1:5
  WrongCellsArrayPrctl{i}=WrongCells;
end;
% allcells = cellfun(@(x,y)x(~ismember(x,y)),allcells,WrongCellsArray,'uni',false);
%%

%  allcells = cellfun(@(x,y)x(~ismember(x,y)),allcells,WrongCellsArray,'uni',false);
% clear all wrong cells out of allcell array 
%%
for i=1:EventsIdxEd-1
    SpEventSize(i)=length(allcells{i+1})/max(cellfun(@length,allcells));
     SpEngCellsRatio(i)=length(intersect(allcells{i+1},allcells{i}))/length(allcells{i+1});
  InterEventDur(i)=EvTime(i+1)-EvTime(i); %time between two consequitive population events
 CommonFraction{i}=intersect(allcells{i+1},allcells{i}); %the common cells participating in both consequitive population events
 CommFrSize(i)=length(CommonFraction{i}); % the number of common cells vectror
 CommFrPercentage(i)=CommFrSize(i)/(length(allcells{i+1})+length(allcells{i})-CommFrSize(i)); %number of common cells weithed by the overall number of cells participating in both events, counting the repeating ones only once, 1 if both events overlap totally, 0 if there is no overlapping at all
 EvSize(i)=length(allcells{i+1})/length(allcells{i}); % the relative size of the following event compared to previous one
 AbsEvSize(i)=length(allcells{i+1})/length(region.onsets);
 EvDurPlot(i)=EvDur(i+1);
end;


EventsIdxEd=EventsIdxEd;
%%
clear('inta')
clear('frb')
j=1;
for i=1:length(CommFrPercentage) 
    if InterEventDur(i)<=100 %choose only for the events with the intervals shorter than 3s or 100 frames
        inta(j)=InterEventDur(i);
        frb(j)= CommFrPercentage(i);
        frbabs(j)=CommFrSize(i);
        evs(j)=EvSize(i);
        Edp(j)=EvDurPlot(i);
        AbEvS(j)=AbsEvSize(i);
        j=j+1;
    end;
end;

[r p]=corrcoef(inta,frb)
save('CommFractionVsIntNoLim','InterEventDur','CommFrPercentage');
save('CommFractionVsInt3sLim','inta','frb');
save('SpVsSpNoIntLim','SpEventSize','SpEngCellsRatio');


%% Plot raster and correlation
clear('EvokedEventIndex')
stim=stim;
fig1=figure;
a1=0.0014;
 b1=-150; %90
%     b1=-3000; %correction for the shifted stimulations
 x(1)=0;
for i=2:length(signals(:,2))
 x(i)=i*T(2).dt+(a1*i*T(2).dt+b1);
% x(i)=i+(0.0016*i);
end; 

StimInd=isempty(stim); %indicator if there is data about the stimulation loaded
if StimInd==0

xstim(1)=0; 
k=1;
StimTimeInd=0;
  StimDuration=500; %for 50 sounds protocol, one stimulus duration is 500 ms
%  StimDuration=150;  %for WN protocol, one stimulus duration is 100 ms

imax=0;
 for i=2:length(stim)
 xstim(i)=i*1; % dt is 1 ms

if (stim(i)>=7)||(stim(i)<=-5) % -5 7,  stim(i)>25 %100 for WN protocol, 25. for 50 sounds protocol 7 -5
    if StimTimeInd==0
    StimTime(k)=i;
    StimTimeInd=1;
    k=k+1;
    imax=i+StimDuration+50; %550 for 50 sounds protocol 150 for WN protocol
    end;
end;
if i>=imax
  StimTimeInd=0;   
end;
 end;  

%   StimTime = [StimTime(1:69) StimTime(69)+1500 StimTime(69)+3000  StimTime(70:end)]; %only one experiment correction, sounds not played
end;
%   
% hax=axes; 


ax(1)=subplot(3,1,1);

 plot(x, Episodes_RealFiltered, 'r')  
 hold on
  plot(x,y,'Color', [0.25, 0.25, 0.25]);
plot(x(EvStart),Episodes_RealFiltered(EvStart),'.','color','bl');
plot(x(EvEnd),Episodes_RealFiltered(EvEnd),'o','color','g');
 if StimInd==0
      plot(xstim,stim/200);
     %===============================================type of stimuli by
     %colour
     
     for k=1:length(StimTime)
        if ((vec.trecord(k)>=1)&&(vec.trecord(k)<=2))||((vec.trecord(k)>=15)&&(vec.trecord(k)<=32))
                       
            if vec.trecord(k)==1||vec.trecord(k)==15||vec.trecord(k)==17||vec.trecord(k)==19||vec.trecord(k)==21||vec.trecord(k)==23||vec.trecord(k)==25||vec.trecord(k)==27||vec.trecord(k)==29||vec.trecord(k)==31
                bar(StimTime(k),max(y)/2,'BarWidth',0.1,'FaceColor','g','EdgeColor','g','LineWidth',0.01)
            else
                bar(StimTime(k),max(y),'BarWidth',0.1,'FaceColor','g','EdgeColor','g','LineWidth',0.01)
            end
        end;
        if ((vec.trecord(k)>=39)&&(vec.trecord(k)<=50))
            if ((vec.trecord(k)>=39)&&(vec.trecord(k)<=44))
             bar(StimTime(k),max(y)/2,'BarWidth',0.1,'FaceColor','r','EdgeColor','r','LineWidth',0.01)
            else
                 bar(StimTime(k),max(y),'BarWidth',0.1,'FaceColor','r','EdgeColor','r','LineWidth',0.01)
            end
        end;
        if ((vec.trecord(k)>=9)&&(vec.trecord(k)<=14))
            if ((vec.trecord(k)>=9)&&(vec.trecord(k)<=11))
             bar(StimTime(k),max(y)/2,'BarWidth',0.1,'FaceColor','k','EdgeColor','k','LineWidth',0.01)
            else
                 bar(StimTime(k),max(y),'BarWidth',0.1,'FaceColor','k','EdgeColor','k','LineWidth',0.01)
            end
        end;
       if ((vec.trecord(k)>=3)&&(vec.trecord(k)<=8))
           if ((vec.trecord(k)>=3)&&(vec.trecord(k)<=5))
             bar(StimTime(k),max(y)/2,'BarWidth',0.1,'FaceColor','m','EdgeColor','m','LineWidth',0.01)
           else
                  bar(StimTime(k),max(y),'BarWidth',0.1,'FaceColor','m','EdgeColor','m','LineWidth',0.01)
           end
        end; 
       if ((vec.trecord(k)>=33)&&(vec.trecord(k)<=38))
             bar(StimTime(k),max(y)/1.5,'BarWidth',0.1,'FaceColor','b','EdgeColor','b','LineWidth',0.01)
        end;   
    end;
     
     %==============================================
 end;    
 plot(x, EventThreshold, 'g')  
%   plot(episodes_res(9,:),'g')

 plot(x,basel,'b')
% plot([0,9543],[mean(episodes_res(96,:)),mean(episodes_res(96,:))]);
%  xlim([0 length(nt)]) 
  xlim([0 max(x)]) 
% plot(episodes_res(5,:),'b')
 
ii=1;
if StimInd==0
 clear('StimCount')
end;
 for i=1:length(EvTime)%length(EventsIdxEd)
%      j=EventsIdxEd(i);
j=i;
    
% line([ min(AboveThEvents{j})  min(AboveThEvents{j})],get(hax,'YLim'),'Color',[0.7 0.7 0.7])
% line([ max(AboveThEvents{j})  max(AboveThEvents{j})],get(hax,'YLim'),'Color',[0.7 0.7 0.7])
txt1 = num2str(i);
% text(min(x(AboveThEvents{j})),1,txt1)
text(x(EvStart(j)),-1,txt1)
if StimInd==0
 
    for k=1:length(StimTime)
        
% if (x(AboveThEvents{j}(1))>=(StimTime(k)-(200-3*(k-1))))&&(x(AboveThEvents{j}(1))<=(StimTime(k)+100+5*k))
if ((x(EvTime(i))>=StimTime(k))&&(x(EvTime(i))<=(StimTime(k)+StimDuration)))||((x(AboveThEvents{j}(1))>=(StimTime(k)-(120-2.5*(k-1))))&&(x(AboveThEvents{j}(1))<=(StimTime(k)+100+5*k)))   %identification of the population events corresponding to auditory response
    EvokedEventIndex(ii)=i; % 500 ms atfer the stimulation to identify the rebound response; based on D. M. Bowman et al 1995 for cata ACx
    DurEvokedEvent(ii)=length(AboveThEvents{j});
    StimCount(ii)=k;
    ii=ii+1;
end
    end;
end;
 end   
%% find the cells responsive to the sound stimulation and sort them according to their responsiveness rate
if StimInd==0
  
   for i=1:length(region.onsets) %scan over each cell onsets
        ii=1;
         StimResponses{i}=[];
          ResponsiveCellRate(i)=0;
            
   
          
     for k=1:length(StimTime) %compare each onsets with stimulation times
         FirstIdx=0;
        for j=1:length(region.onsets{i}) 
            if FirstIdx==0;
           if ((x(region.onsets{i}(j))>=StimTime(k))&&(x(region.onsets{i}(j))<=(StimTime(k)+StimDuration))) %% very important stimdulation +500, 22 avr 2020
               ResponsiveCellRate(i)=ii;
               StimResponses{i}(ii)=k;
               ii=ii+1;
               FirstIdx=1;
           end;
            end
        end;
        
     end;

   end;
    
   
   
   
    [RC,I] = sort(ResponsiveCellRate,'descend');
    % Range of neurons by their responsiveness
ASR=cell2mat(StimResponses); %vector array from allcells containing the neurons participating in each population event
T = tabulate(ASR); %attribute to each neuron its participation rate
[~,idx] = sort(T(:,2),'descend'); % sort the neurons according to their participation rate, from min to max
sortedStim = T(idx,:); % sort thestimulation according to the number of neurons responding
sortedStim(1); %the most evoking stimulation
%% neurons responding to each stimulation
 for i=1:length(StimTime)
     ii=1;
     NeuronsPerStim{i}=[];
  for  k=1:length(StimResponses)
      if isempty(StimResponses{k})==0
      for j=1:length(StimResponses{k})
      if StimResponses{k}(j)==i
NeuronsPerStim{i}(ii)=k;    %cell arrays of neurons responding to each auditory stimulation
ii=ii+1;
      end;
      end;
      end;
  end;
  NeuronsPerStim{i}=unique(NeuronsPerStim{i});
 end;
    
end;
%%
   hold off
  
   title('Population firing rate and synchronous events')


ax(2)=subplot(3,1,2);
 palette=distinguishable_colors(length(region.onsets));
% RAW TRACES PLOT
%   plot(Episodes_RealFiltered, 'r')  
%  plot(Episodes_RealFiltered, 'r')  
%   plot(nt(1,:))
   hold on
  step=0;
  
  for i=1:134 %length(ParticipativeCells0)%length(region.onsets) Plotting the traces for the 3 most participative cells
      step=step+100;
      plot(x(:),nt(sortedmat(i),:),'color',palette(i,:,:))


  end;
      xlim([0 max(x)]) 
 
   hold off
% STIM or EYE Pupil traces plot
%      plot(stim)
EyeInd=exist('allres','var'); 
if EyeInd==1
    title('pupil tracking')
         plot(allres.pupil(:,3))
         xlim([0 length(allres.pupil(:,3))]) 
end;
%            xlim([0 length(stim)]) 
%      title('white noise stimulation')


%  title('raw traces')

%  fig2=figure;

ax(3)=subplot(3,1,3);
for itrans=1:length(region.onsets)
    for jtrans=1:length(region.onsets{itrans})
    regionMod.onsets{itrans}(jtrans)=x(region.onsets{itrans}(jtrans));
    end;
    if length(region.onsets{itrans})==0
      regionMod.onsets{itrans}=[0];       
    end;
end;    
plotSpikeRaster(regionMod.onsets);
title('Raster plot')
% imagesc(all_s{1});
linkaxes(ax,'x'); 
fig2=figure;
for i=1:length(allcells)
    EventSize(i)=length(allcells{i});
    if StimInd==0
   for j=1:length(EvokedEventIndex)
       if i==EvokedEventIndex(j)
           EvokedEvSize(j)=length(allcells{i});
          
       end
   end;
    end;  
    
    
end;
% hist(EventSize,round(max(EventSize)/10))
 hist(AbsEvSize)
 if StimInd==0
hist(EvokedEvSize/length(region.onsets))
end;
%% 
if StimInd==0
    for i=1:length(EvokedEventIndex)
 EvokedEvDur(i)=x(EvEnd(EvokedEventIndex(i)))- x(EvStart(EvokedEventIndex(i)));
     EvWindowStart(i)=EvStart(EvokedEventIndex(i))-30;
     EvWindowEnd(i)=EvStart(EvokedEventIndex(i));
    end;
    end;
    
 for   i=1:length(EvTime) 
     SpEvDur(i)=x(EvEnd(i))-x(EvStart(i));
     WindowStart(i)=EvStart(i)-30;
     WindowEnd(i)=EvStart(i);
 end;
 
%%
for i=1:size(region.traces,1)
   FR(i)=length(region.onsets{i})/300;
    CV(i)=std(diff(region.onsets{i}))/mean(diff(region.onsets{i}));
    Fano(i)=var(diff(region.onsets{i}))/mean(diff(region.onsets{i}));
%     diff(region.onsets{i});
end;    
nanmean(CV) %coefficient of variation of the interspike interval
nanmean(FR) %average firing rate
nanmean(Fano) %average Fano factor


%% Plot of the participative cells in all population events
% BB=length(sortedmat(:,1,1));
% CellCountAll=ParticipationThresholdSp(numres,ArtSpTr,EvTime,EvEnd,EvStart);
palette=jet(sortedmat(1,2,1));
%  B = fliplr(palette);
% palette1=sort(palette,'descend');
fig3=figure;
 title('The most participative cells >50% of population events')
 ii=1;
for i=1:length(sortedmat(:,1,1))
    hold on
    if sortedmat(i,2,1)>mean(CellCount); %(sortedmat(1,2,1)/2)
 hippo_plot_cont_clustering(region, sortedmat(i,1,1), palette(sortedmat(i,2,1),:));
 ParticipativeCells0(ii)=sortedmat(i,1,1);
 ParticipationRate0(ii)=sortedmat(i,2,1);
  ii=ii+1;
    end
end;
 
 ii=1;
for i=1:length(ParticipativeCells0)
   
    if ParticipationRate0(i)>prctile(ParticipationRate0,95) %(sortedmatSp(1,2,1)/2)
       ParticipativeCells095(ii)=ParticipativeCells0(i);
       ii=ii+1;
    end;
end;

ParticipativeCells095=setdiff(ParticipativeCells095,WrongCells);
NumNeurons=size(region.traces,1);
save([fullfile(pwd,'\') 'SpontOnly' num2str(sscanf(filename,'region%d')) '.mat'],'EvTime','allcells','NumNeurons','EvStart','EvEnd');
% save('SpontOnly','EvTime','allcells','NumNeurons','EvStart','EvEnd');
% fig6=figure;
% for i=1:3
%     hold on
%  hippo_plot_cont_clustering(region, sortedmat(i,1,1), B(i,:));
% end;
%% 
if StimInd==0 %in case there is any audotory stimulation data
%% Only evoked events
Aev=cell2mat(allcells(EvokedEventIndex)); %vector array from allcells containing the neurons participating in each evoked event
Tev = tabulate(Aev); %attribute to each neuron its participation rate
Tev=Tev(find(Tev(:,2,1)>=1),:,:);
[~,idx] = sort(Tev(:,2),'descend'); % sort the neurons according to their participation rate, from min to max
sortedmatEv = Tev(idx,:); % sort the neurons according to their participation rate, from min to max
sortedmatEv(1); %the most participative neuron
% BBev=length(sortedmatEv(:,1,1));
palette=jet(sortedmatEv(1,2,1));
% palette=jet(BBev);
B = fliplr(palette);
fig7=figure;
for i=1:length(sortedmatEv(:,1,1))
    hold on
 hippo_plot_cont_clustering(region, sortedmatEv(i,1,1), palette(sortedmatEv(i,2,1),:));
end;

fig8=figure;
for i=1:20
    hold on
 hippo_plot_cont_clustering(region, sortedmatEv(i,1,1), palette(sortedmatEv(i,2,1),:));
end;
%% Only spontaneous event

startValue = 1;
endValue = length(allcells);
nElements = length(allcells);
A = linspace(startValue,endValue,nElements);

OnlySpontInd=setdiff(A,EvokedEventIndex);
ASp=cell2mat(allcells(OnlySpontInd)); %vector array from allcells containing the neurons participating in each only spontaneous population event
TSp = tabulate(ASp);%attribute to each neuron its participation rate
TSp=TSp(find(TSp(:,2,1)>=1),:,:);
[~,idx] = sort(TSp(:,2),'descend'); % sort the neurons according to their participation rate, from min to max
sortedmatSp = TSp(idx,:); % sort the neurons according to their participation rate, from min to max
% sorting of the cells according to their lattency from the beginning of
% the spontaneous event

   CellCountSp=ParticipationThresholdSp(numres,ArtSpTr(sortedmatSp(find(sortedmatSp(:,2)>0),1)),EvTime(EvokedEventIndex),EvEnd(EvokedEventIndex),EvStart(EvokedEventIndex));

palette=jet(sortedmatSp(1,2,1));
%  B = fliplr(palette);
% palette1=sort(palette,'descend');
fig9=figure;
 title('The most participative cells >50% of only spontaneous evens')
  ParticipativeCellsPrctile=cell(1,5);
 ii=1;
for i=1:length(sortedmatSp(:,1,1))
    hold on
    if sortedmatSp(i,2,1)>mean(CellCountSp); %(sortedmatSp(1,2,1)/2)
 hippo_plot_cont_clustering(region, sortedmatSp(i,1,1), palette(sortedmatSp(i,2,1),:));
  ParticipativeCells(ii)=sortedmatSp(i,1,1);
  ParticipationRate(ii)=sortedmatSp(i,2,1);
  ii=ii+1;
    end
end;
%%  
%   ParticipativeCellsPrctile{1}=ParticipativeCells;
% ii=1;
% ii2=1;
% ii3=1;
% ii4=1;
% ii5=1;
% for i=1:length(ParticipativeCells)
%     if ParticipationRate(i)>prctile(ParticipationRate,25) %(sortedmatSp(1,2,1)/2)
%        ParticipativeCellsPrctile{2}(ii2)=ParticipativeCells(i);
%        ii2=ii2+1;
%     end;
%       if ParticipationRate(i)>prctile(ParticipationRate,50) %(sortedmatSp(1,2,1)/2)
%        ParticipativeCellsPrctile{3}(ii3)=ParticipativeCells(i);
%        ii3=ii3+1;
%      end;  
%    if ParticipationRate(i)>prctile(ParticipationRate,75) %(sortedmatSp(1,2,1)/2)
%        ParticipativeCellsPrctile{4}(ii4)=ParticipativeCells(i);
%        ii4=ii4+1;
%      end;    
%     if ParticipationRate(i)>prctile(ParticipationRate,95) %(sortedmatSp(1,2,1)/2)
%        ParticipativeCellsPrctile{5}(ii5)=ParticipativeCells(i);
%        ii5=ii5+1;
%      end;      
% end;

% ParticipativeCells95=setdiff(ParticipativeCells95,WrongCells);

% ParticipativeCellsPrctile = cellfun(@(x,y)x(~ismember(x,y)),ParticipativeCellsPrctile,WrongCellsArrayPrctl,'uni',false);
% 
% fig10=figure;
% for i=1:20
%     hold on
%  hippo_plot_cont_clustering(region, sortedmatSp(i,1,1), B(i,:));
% end;
%% 
MeanSizeOnlySpEvent = mean(cellfun(@length, allcells(OnlySpontInd)));
MeanSizeOnlyEvokedEvent = mean(cellfun(@length, allcells(EvokedEventIndex)));
end;

%% check of activity before the evoked events
PrevWindow=1000; %ms
PrWinSec=1; %in seconds
if StimInd==0 %in case there is any audotory stimulation data
    for i=1:length(StimTime)
         PrevEvCount(i)=0;
        for j=1:length(EvTime)
            if (x(EvTime(j))<StimTime(i))&&(x(EvTime(j))>(StimTime(i)-PrevWindow)) %if there is any population event within 1 second before the stimulation
                %nothing, one or more than one population event
                PrevEvCount(i)=PrevEvCount(i)+1;
            end;
        end;
    end;
  
    
for i=1:length(StimTime)
    FRPrevCell=[];
    PrevSpks=[];
    for j=1:length(all_s{1}(:,1))
        PrevSpks=x(all_s{1}(j,:)==1);
        FRPrevCell(j)=length( find((PrevSpks>(StimTime(i)-PrevWindow))&(PrevSpks<=StimTime(i))))/PrWinSec;
    end;
   
    FRprevPop(i)=mean(FRPrevCell); % average FR for individual cells during 3s before the stimulation
end;
   
   for i=1:length(StimTime)
       StimIdx=find(StimCount==i);
       if isempty(StimIdx)==1
           RespAmpl(i)=0;
       else
     RespAmpl(i)=length(allcells{EvokedEventIndex(StimIdx(1))});
       end;

   end;
    
end;
%%  Responsive cells plot
clear('ResponsiveCells')
if StimInd==0
       CellCountEvoked=ParticipationThresholdEvoked(x,numres,ArtSpTr(sortedmatEv(find(sortedmatEv(:,2)>0),1)),StimTime,StimTime+StimDuration,StimTime);
palette=jet(RC(1));
% B = fliplr(palette);
fig10=figure;

   title('The most responsive cells >50% of stimulations')
  ResponsiveCellsPrctile=cell(1,5);
   ii=1;
for i=1:length(I)
      hold on
    if RC(i)>mean(CellCountEvoked); %RC(1)/2 %the most responsive cells, responding to more than a half of the highest rate. 
      
 hippo_plot_cont_clustering(region, I(i), palette(round(RC(i)),:));
 ResponsiveCells(ii)=I(i);
 ResponsiveRate(ii)=RC(i);
 ii=ii+1;
    end;
end;
%% 
%  ResponsiveCellsPrctile{1}=ResponsiveCells;
% ii=1;
% ii2=1;
% ii3=1;
% ii4=1;
% ii5=1;
% for i=1:length(ResponsiveCells)
%      
%     if ResponsiveRate(i)>prctile(ResponsiveRate,25); %RC(1)/2 %the most responsive cells, responding to more than a half of the highest rate. 
%               ResponsiveCellsPrctile{2}(ii2)= ResponsiveCells(i);
%          %  ResponsiveRate(ii)=RC(i);
%               ii2=ii2+1;
%     end;
%   if ResponsiveRate(i)>prctile(ResponsiveRate,50); %RC(1)/2 %the most responsive cells, responding to more than a half of the highest rate. 
%       ResponsiveCellsPrctile{3}(ii3)= ResponsiveCells(i);
%       %  ResponsiveRate(ii)=RC(i);
%       ii3=ii3+1;
%   end;
%      if ResponsiveRate(i)>prctile(ResponsiveRate,75); %RC(1)/2 %the most responsive cells, responding to more than a half of the highest rate. 
%        ResponsiveCellsPrctile{4}(ii4)= ResponsiveCells(i);
%        %  ResponsiveRate(ii)=RC(i);
%        ii4=ii4+1;
%      end;
%      if ResponsiveRate(i)>prctile(ResponsiveRate,95); %RC(1)/2 %the most responsive cells, responding to more than a half of the highest rate. 
%       ResponsiveCellsPrctile{5}(ii5)= ResponsiveCells(i);
%       %  ResponsiveRate(ii)=RC(i);
%       ii5=ii5+1;
%     end;
% end;
% 
% ResponsiveCellsPrctile = cellfun(@(x,y)x(~ismember(x,y)),ResponsiveCellsPrctile,WrongCellsArrayPrctl,'uni',false);
%  ResponsiveCells95=setdiff(ResponsiveCellsPrctile,WrongCells);
%% 
% length(intersect(ResponsiveCells95,ParticipativeCells95))/(length(ResponsiveCells95)+length(ParticipativeCells95)-length(intersect(ResponsiveCells95,ParticipativeCells95)))
% for i=1:5
% OverlapEvSp(i)=length(intersect(ResponsiveCellsPrctile{i},ParticipativeCellsPrctile{i}))/(length(ResponsiveCellsPrctile{i})+length(ParticipativeCellsPrctile{i})-length(intersect(ResponsiveCellsPrctile{i},ParticipativeCellsPrctile{i})))
% end;
NumNeurons=size(region.traces,1);
 SoundIndex=vec.trecord;
 save([fullfile(pwd,'\') 'EvokedVsSpontStim' num2str(sscanf(filename,'region%d')) '.mat'],'EvTime','StimTime','OnlySpontInd','NeuronsPerStim','allcells','NumNeurons','SoundIndex');
% save('EvokedVsSpontStim','EvTime','StimTime','OnlySpontInd','NeuronsPerStim','allcells','NumNeurons','SoundIndex');
end;


 
%% 
%% 
% fig5=figure;
%  level=0;
%  for i=1:length(all_cluster)
%   %   colours(i,1:3)=palette(round(BB(i,1)/(j-1)*100),1:3);
%         hold on 
%         if length(all_cluster{:,i})>1
%         rasterpoint(all_cluster,region1.onsets(all_cluster{i}),level,palette(i,:))
%         level=level+length(all_cluster{i});
%         end;
%  end;