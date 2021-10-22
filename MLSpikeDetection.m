NumOfFrames=length(signals(:,1)); 
p = spk_est('par');
for i=1:length(regioncenters)
    





% p = struct('output','MAP spike train', ...
%                 'a',0.1,'tau',0.76,'ton',0, ...
%                 'nonlinearity',[], ...
%                 'saturation',0.1,'hill',2.99, ...
%                 'p2',0.5,'p3',0.01, ...
%                 'sigma',[],'drift',.015, ...
%                 'bmin',[],'bmax',[], ...
%                 'cmax',10,'nc',50,'np',50,'nb',50);
% p.a=str2num(get(a,'string'));
p.dt=300/NumOfFrames;
p.tau=str2num(get(tau,'string'));
p.pnonlin=[];
% p.saturation = str2num(get(saturation,'string'));
% p.hill=str2num(get(hill,'string'));
p.pnonlin = [.5 .01];
p.drift.parameter=str2num(get(drift0,'string')); %0.015;
p.finetune.sigma = []; % auto estimation of sigma, this is already the default value
calcium=transpose(signals(:,i,TrNum)); 
dt=p.dt; 





% positionVector1 = [0.03, 0.50, 0.94, 0.17];    % position of first subplot
% subplot('Position',positionVector1)
%   plot(calcium); 
%   axis([0 9543 0 inf])
%   hold on
%   hold off

j=i;
positionVector3 = [0.03, 0.06, 0.94, 0.17];    % position of second subplot
subplot('Position',positionVector3)
[spikest fit drift] = spk_est(calcium,p);

              Y(1:length(spikest))=min(calcium);%i*2; 
             
              plot(fit,'r');  
              axis([0 NumOfFrames min(calcium)-5 inf])
              hold on
             scatter(spikest/dt,Y,10,'filled');
             
            hold off    
         
              Y=[];
              i
              

set(txcellnum, 'String', num2str(i)); 
positionVector1 = [0.03, 0.50, 0.94, 0.17];    % position of first subplot
subplot('Position',positionVector1)
plot(calcium); 
axis([0 NumOfFrames 0 inf])

 hold off

%%
num=i;
eq=0;
 f = find(spk(num,:)==1);
 g = find(dec(num,:)==1);


if numOld~=num
    subplot(imgax)
    hold off
%     imagesc(a)
%     hold on
    set(gca,'xtick',[],'ytick',[]);
    axis equal
    axis tight
    box on
    colormap gray
    set(bbright,'enable','on');
    set(bcontrast,'enable','on');
    [maxy maxx] = size(a);
    HippoContrast;
    cellSel=num;
    for i=cellSel
        hold on;
          plot(region.contours{i}([1:end 1],1),region.contours{i}([1:end 1],2),'r-')
    end
end
subplot(trax)
%% 
 
 region.onsets{j}=round(spikest/dt);
 region.offsets{j}=round(spikest/dt);    
 region.fits(j,:)=fit;
end;   

region.tau=str2num(get(tau,'string'));
region.drift=str2num(get(drift0,'string'));
% save(['/mnt/storage05EqBathellier_2w/ANALYSIS/Anton/200415/region ' num2str(TrNum) '.mat'],'region');
save([SaveRegion num2str(TrNum) '.mat'],'region');
clear
     save()         
        