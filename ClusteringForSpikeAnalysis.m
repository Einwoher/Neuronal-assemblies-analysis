%    big_cluster1=[allcells1 allcells2 allcells3];
%  big_cluster1=[allcells(OnlySpontInd)];
big_cluster1=[allcells];
allcells=[];
 allcells= big_cluster1;
%  
% soundindex=[SoundIndex1,SoundIndex2,SoundIndex3];
% NeuronsPerStim=[NeuronsPerStim1,NeuronsPerStim2,NeuronsPerStim3];
% allcells=NeuronsPerStim;
% big_cluster1=NeuronsPerStim;

C=big_cluster1;
for i=1:length(C) 
    for j=1:length(C)
%         a(i,j)=[];
       




%         a(i,j)=length(intersect(cell2mat(C(i)),cell2mat(C(j))))/length(union(cell2mat(C(i)),cell2mat(C(j)))); %Jaccard similarity
%          a(i,j)=length(intersect(cell2mat(C(i)),cell2mat(C(j))))/ min(length(cell2mat(C(i))),length(cell2mat(C(j)))); % Szymkiewicz–Simpson coefficient
%            a(i,j)=2*length(intersect(cell2mat(C(i)),cell2mat(C(j))))/(length(cell2mat(C(i)))+length(cell2mat(C(j)))); %Sorensen Coefficient
vec1=[];
vec2=[];
  vec1=zeros(NumNeurons,1); % model
  vec2=zeros(NumNeurons,1); % model
%  vec1=zeros(length(region.onsets),1);
%  vec2=zeros(length(region.onsets),1);
 
  positions1=C{i};
  positions2=C{j}; 
 vec1(positions1)=1;
 vec2(positions2)=1;
 R=corrcoef(vec1,vec2);  
  a(i,j)=R(1,2); %Pearsons coefficient

%         (length(C{i})+length(C{j})-length(intersect(cell2mat(C(i)),cell2mat(C(j)))));
           
    end;
end;
% a(isnan(a))=0;
%%
%  a=c;
figure
imagesc(a)
hold on
line([0 length(C)+0.5],[length(big_cluster1)+0.5 length(big_cluster1)+0.5],'Color','k')
% line([0 length(C)+0.5],[length(big_cluster1)+length(big_cluster2)+0.5 length(big_cluster1)+length(big_cluster2)+0.5],'Color','k')
% line([0 length(C)+0.5],[length(big_cluster1)+length(big_cluster2)+length(big_cluster3)+0.5 length(big_cluster1)+length(big_cluster2)+length(big_cluster3)+0.5],'Color','k')
% line([0 length(C)+0.5],[length(big_cluster1)+length(big_cluster2)+length(big_cluster3)+length(big_cluster4)+0.5 length(big_cluster1)+length(big_cluster2)+length(big_cluster3)+length(big_cluster4)+0.5],'Color','k')
% line([0 length(C)+0.5],[length(big_cluster1)+length(big_cluster2)+length(big_cluster3)+length(big_cluster4)+length(big_cluster5)+0.5 length(big_cluster1)+length(big_cluster2)+length(big_cluster3)+length(big_cluster4)+length(big_cluster5)+0.5],'Color','k')
% line([0 length(C)+0.5],[length(big_cluster1)+length(big_cluster2)+length(big_cluster3)+length(big_cluster4)+length(big_cluster5)+length(big_cluster6)+0.5 length(big_cluster1)+length(big_cluster2)+length(big_cluster3)+length(big_cluster4)+length(big_cluster5)+length(big_cluster6)+0.5])

line([length(big_cluster1)+0.5 length(big_cluster1)+0.5], [0 length(C)+0.5],'Color','k')
% line([length(big_cluster1)+length(big_cluster2)+0.5 length(big_cluster1)+length(big_cluster2)+0.5],[0 length(C)+0.5],'Color','k')
% line([length(big_cluster1)+length(big_cluster2)+length(big_cluster3)+0.5 length(big_cluster1)+length(big_cluster2)+length(big_cluster3)+0.5],[0 length(C)+0.5],'Color','k')
% line([length(big_cluster1)+length(big_cluster2)+length(big_cluster3)+length(big_cluster4)+0.5 length(big_cluster1)+length(big_cluster2)+length(big_cluster3)+length(big_cluster4)+0.5],[0 length(C)+0.5],'Color','k')
% line([length(big_cluster1)+length(big_cluster2)+length(big_cluster3)+length(big_cluster4)+length(big_cluster5)+0.5 length(big_cluster1)+length(big_cluster2)+length(big_cluster3)+length(big_cluster4)+length(big_cluster5)+0.5],[0 length(C)+0.5],'Color','k')
% line([length(big_cluster1)+length(big_cluster2)+length(big_cluster3)+length(big_cluster4)+length(big_cluster5)+length(big_cluster6)+0.5 length(big_cluster1)+length(big_cluster2)+length(big_cluster3)+length(big_cluster4)+length(big_cluster5)+length(big_cluster6)+0.5],[0 length(C)+0.5])

%% 

%# remove diagonal elements
corrMat=a;
corrMat = corrMat - eye(size(corrMat));
%# and convert to a vector (as pdist)
% dissimilarity = 1 - corrMat(find(corrMat))';
%# decide on a cutoff
%# remember that 0.4 corresponds to corr of 0.6!
  cutoff =1; 


%# perform complete linkage clustering
% Z = linkage(dissimilarity,'complete');

Z = linkage(corrMat,'complete','correlation');
I = inconsistent(Z);
% cutoff=median(I(:,4)); % medial of inconsistency index used as a cutoff threshold
%# group the data into clusters
%# (cutoff is at a correlation of 0.5)
groups = cluster(Z,'cutoff',cutoff,'criterion','distance');
% figure
% dendrogram(Z)  
u=unique(groups);
for i=1:length(u)
    gr{i}=[];
for j=1:length(a)
   if groups(j)==u(i)
       gr{i}=[gr{i},j];
   end; 
end;
end;
b=cell2mat(gr);

% plot clustered correlation matrix
% C=[];
% C=allcells(b);
% for i=1:length(C) 
%     for j=1:length(C)
% %         a(i,j)=[];
%        
% 
% 
% 
% 
% %        a(i,j)=length(intersect(cell2mat(C(i)),cell2mat(C(j))))/length(union(cell2mat(C(i)),cell2mat(C(j)))); %Jaccard similarity
% %          a(i,j)=length(intersect(cell2mat(C(i)),cell2mat(C(j))))/ min(length(cell2mat(C(i))),length(cell2mat(C(j)))); % Szymkiewicz–Simpson coefficient
% %            a(i,j)=2*length(intersect(cell2mat(C(i)),cell2mat(C(j))))/(length(cell2mat(C(i)))+length(cell2mat(C(j)))); %Sorensen Coefficient
% vec1=[];
% vec2=[];
%   vec1=zeros(NumNeurons,1); % model
%   vec2=zeros(NumNeurons,1); % model
% %  vec1=zeros(length(region.onsets),1);
% %  vec2=zeros(length(region.onsets),1);
%  
%   positions1=C{i};
%   positions2=C{j}; 
%  vec1(positions1)=1;
%  vec2(positions2)=1;
%  R=corrcoef(vec1,vec2);  
%   a(i,j)=R(1,2); %Pearsons coefficient
% 
% %         (length(C{i})+length(C{j})-length(intersect(cell2mat(C(i)),cell2mat(C(j)))));
%            
%     end;
% end;

figure
imagesc(a(b,b))
ClusterColour=distinguishable_colors(length(gr));
figure
imagesc(a(b,b))
hold on
 Position=0;
 for i=1:length(gr)
   rectangle('Position',[Position,Position,length(gr{i}), length(gr{i})],'Edgecolor', ClusterColour(i,:,:),'FaceColor',ClusterColour(i,:,:));
  text('Position', [Position,Position], 'String', num2str(i));
  num2str(i)
  Position=Position+length(gr{i});
 end 
 save([fullfile(pwd,'\') 'EventClusters' num2str(sscanf(filename,'region%d')) '.mat'],'b','gr','allcells','NumNeurons');
%    save('EventClusters2','b','gr','allcells','soundindex');
%      save('EventClusters','b','gr','allcells','NumNeurons');