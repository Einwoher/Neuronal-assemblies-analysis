
function NeuropilCorr = NeuropilCorr(signals,localneuropil)
signals1=signals-0.7*localneuropil;

% Correction for feeble signal
% newsignal0=signals1*3;
% CorCoef=abs(mean(newsignal0)-mean(signals1))
% newsignal=newsignal0-CorCoef;
% signals1=newsignal0;
%end of correction

a=min(signals1);
  b=mean(0.7*localneuropil);
for i=1:length(signals1(1,:))
    signals1(:,i)=signals1(:,i)+abs(a(i))+abs(b(i));
%      signals1(:,i)=signals1(:,i)+abs(a(i));
end;
%  signals=signals1;
NeuropilCorr= signals1;
% signals2=signals';
% 
% 
% for c = 1:size(signals2,1)
%      signals3(c,:) = dfoverf(signals2(c,:))*100;
%     
% end;
% 
% signals=signals3';
%  a=min(signals);
% % 
% for i=1:length(signals(1,:))
%     signals(:,i)=signals(:,i)+abs(a(i));
% end;