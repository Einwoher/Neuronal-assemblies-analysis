function dfoverf = dfoverf(tr)

       bases = repmat(mean(tr,2),1,size(tr,2)); %old verion with just a
%      average baseline
  y=tr;
%  lambda=1.0000e+09;
%      p= 0.0100;
% 
% y=double(tr);
%     bases=baseline_(y',lambda,p);
% 
%   bases = bases';
dfoverf = (y-bases)./(bases+eps);



% if mean(f) == 0
%    dfoverf = f / (max(f)-min(f)+eps);
% else
%    a = (f - median(f))/median(f);
%    dfoverf = a ;
% end
