%%
if strcmp(get(EventsToPlot,'string'),'all')==1
  EndEvent=length(allcells);
  StartEvent=1;
elseif str2double(get(EventsToPlot,'string'))>(length(allcells))
    set(EventsToPlot, 'String', num2str(length(allcells))); 
    EndEvent=length(allcells);
    StartEvent=length(allcells);
else
   EndEvent=str2double(get(EventsToPlot,'string'));
   StartEvent= EndEvent;
end
    
  fig4=figure;
for im=StartEvent:EndEvent
EvNum=0+im; %number of the event to plot
EventsIdx(EvNum)
region1=region;
 txt2 = num2str(im);
% for i = 1:size(region1.traces,1)
%     all_cells{1}(i,region1.onsets{i})=1;% creates a binary matrix with 0 and 1, where each event onset is 1
%     all_cells{2}(i)=i; 
% end
%  take off the waves

lim1=min(AboveThEvents{EventsIdx(EvNum)});
lim2=max(AboveThEvents{EventsIdx(EvNum)});





% Wave propagation
 clear('all_cluster')
 clear('palette');
 all_cluster={};
startwave=lim1;
endwave=lim2;
for i=startwave:endwave
     k=0;
    for j=1:length(region1.onsets)
      if ismember(i,region1.onsets{j})==1
        k=k+1;
        all_cluster{i-startwave+1}(k)=j;
      end;
    end;
  end;

BB=length(all_cluster);
palette=jet(BB);
% palette(1,:,:)=[1 0 0];
% palette(length(all_cluster),:,:)=[0 0 0.55];
 palette(1,:,:)=[0 0 0.55];
 palette(length(all_cluster),:,:)=[1 0 0];
palette=flipud(palette); %start with hot colours




for i = 1:size(region1.traces,1)
     for j=1:length(region1.onsets{i});
lim1=min(AboveThEvents{EventsIdx(EvNum)});
lim2=max(AboveThEvents{EventsIdx(EvNum)});

if (region1.onsets{i}(j)<lim1)||(region1.onsets{i}(j)>lim2) %wave1 condtion
         region1.onsets{i}(j)=0;
         region1.offsets{i}(j)=0;
      
end    


     end
end
% for i=1:size(region1.traces,1)
%    a=region1.onsets{i};
%    b=region1.offsets{i};
%    jj=length(a);
%    j=1;
%    while j<=jj
%    if a(j)==0 
%        a(j)=[];
%        jj=jj-1;
%        j=0;
%    end
%    j=j+1;
%    end
%    jj=length(b);
%    j1=1;
%    while j1<=jj
%    if b(j1)==0 
%        b(j1)=[];
%        jj=jj-1;
%        j1=0;
%    end
%    j1=j1+1;
%    end
%    
%    region1.onsets{i}=a;
%     region1.offsets{i}=b;
%   
% end

 
 
 
% subplot(1,1,1)
if strcmp(get(EventsToPlot,'string'),'all')==1
% subplot(10,16,im)
% FigEvent=figure;
figure;
  subplot(1,1,1);
else
    subplot(1,1,1)
end;
          title(txt2,'FontSize',5)

 for i=1:length(all_cluster)
  %   colours(i,1:3)=palette(round(BB(i,1)/(j-1)*100),1:3);
        hold on
        if length(all_cluster{:,i})>1   
                 hippo_plot_cont_clustering(region1, all_cluster{i}, palette(i,:));
        end;
    
 end;
%  saveas(fig,'txt2.png')
set(gcf, 'InvertHardCopy', 'off');
% fig.PaperPositionMode = 'manual';
%   saveas(gcf,['event' txt2 '.png']);
%  
 



%  if strcmp(get(EventsToPlot,'string'),'all')==0
%  fig5=figure;
%  level=0;
%  for i=1:length(all_cluster)
%   %   colours(i,1:3)=palette(round(BB(i,1)/(j-1)*100),1:3);
%         hold on 
%         if length(all_cluster{:,i})>1
%         rasterpoint(all_cluster{i},region1.onsets(all_cluster{i}),level,palette(i,:))
%         level=level+length(all_cluster{i});
%         end;       
%  end;
%  end;
end;

