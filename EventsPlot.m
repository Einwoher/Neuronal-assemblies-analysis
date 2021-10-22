% niter=100; sig=1; p_val=0.95; %niter number of iterations for surrogate data set
all_cells = cell(1,2); % creates a cell array
all_cells{1} = zeros(size(region.traces)); % all neuron traces 
all_cells{2} = zeros(1,size(region.traces,1)); % corresponding cell numbers
lim1ms=str2double(get(EventStart,'string')); %from ms to frames  
lim1=round((lim1ms-b1)/(T(2).dt*(1+a1)));

lim2ms=str2double(get(EventEnd,'string')); %%from ms to frames
lim2=round((lim2ms-b1)/(T(2).dt*(1+a1)));
for i = 1:size(region.traces,1)
    all_cells{1}(i,region.onsets{i})=1;% creates a binary matrix with 0 and 1, where each event onset is 1
    all_cells{2}(i)=i; 
end
region.onsetsB=region.onsets;
%  take off the waves
for i = 1:size(region.traces,1)
     for j=1:length(region.onsets{i});

if (region.onsetsB{i}(j)<lim1)||(region.onsetsB{i}(j)>lim2) %wave1 condtion
         region.onsetsB{i}(j)=0;
         region.offsetsB{i}(j)=0;
end    


     end
end
for i=1:size(region.traces,1)
   a=region.onsetsB{i};
   b=region.offsetsB{i};
   jj=length(a);
   j=1;
   while j<=jj
   if a(j)==0 
       a(j)=[];
       jj=jj-1;
       j=0;
   end
   j=j+1;
   end
   jj=length(b);
   j1=1;
   while j1<=jj
   if b(j1)==0 
       b(j1)=[];
       jj=jj-1;
       j1=0;
   end
   j1=j1+1;
   end
   
   region.onsetsB{i}=a;
    region.offsetsB{i}=b;
  
end
% Wave propagation
clear('all_cluster')
all_cluster={}
startwave=lim1;
endwave=lim2;
for i=startwave:endwave
     k=0;
    for j=1:length(region.onsetsB)
      if region.onsetsB{j}==i
        k=k+1;
        all_cluster{i-startwave+1}(k)=j;
      end;
    end;
end;
BB=length(all_cluster);
palette=jet(BB);
 palette=flipud(palette); %start with hot colours
 
%% 
 fig1=figure;
 subplot(1,2,1) 
 for i=1:length(all_cluster)
  %   colours(i,1:3)=palette(round(BB(i,1)/(j-1)*100),1:3);
        hold on
        if length(all_cluster{:,i})>1          
        hippo_plot_cont_clustering(region, all_cluster{i}, palette(i,:))
        end;
    
 end;
title('Spatial propagation map')
subplot(1,2,2) 
 level=0;
 for i=1:length(all_cluster)
  %   colours(i,1:3)=palette(round(BB(i,1)/(j-1)*100),1:3);
        hold on 
        if length(all_cluster{:,i})>1
        rasterpoint(all_cluster,region.onsetsB(all_cluster{i}),level,palette(i,:))
        level=level+length(all_cluster{i});
        end;
 end;
 title('Avalanche temporal propagation')