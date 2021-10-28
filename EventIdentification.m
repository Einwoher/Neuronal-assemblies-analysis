%EVENT IDENTIFICATION FUNCTION
% clear('region.onsets');
% clear('region.offsets');
onsetmat(num,1:length(nt(1,:)))=0;
offsetmat(num,1:length(nt(1,:)))=0;
spk(num,1:length(nt(1,:)))=0;
dec(num,1:length(nt(1,:)))=0;

%for num=1:size(nt,1)
%baseline plot:
%extrema calculation
clear('ymax');
clear('ymin');
clear('imin');
clear('imax');
[ymax,imax,ymin,imin] = extrema(nt(num,:)); % num is the number of the actual cell cell N
hold on
%plot(imax,ymax,'r*',imin,ymin,'g*') %extrema plot
%plot(imax,ymax,'r',imin,ymin,'g') %extrema plot
% t = imax(1):200:size(region.traces,2)+imax(1)+1;
% p = pchip(imax,ymax,t);
%nt0=nt(num,:);
%ycorr = bf(nt0,5); 
 %ycorr = msbackadj(1:999,nt(num,:)); 
% plot (ycorr,'r');
%plot (t,p,'r');
% t = imin(1):200:size(region.traces,2)+imin(1)+1;
% p1 = pchip(imin,ymin,t);

%plot (t,p1,'g');

%  noise = mean(abs(p-p1));
%  line = p1+(abs(p-p1))/2;
% plot(t, line+noise,'r');
% plot(t, line-noise,'g');

% P=findvalleys_(1:999,nt(num,:),0.0003,5,13,11,3);
% if P(:,2)>0
% plot(floor(P(:,2))',nt(num,floor(P(:,2))'),'r*');
% end;

thhold=get(threshold,'value');
thhold=thhold*10;
eventdur=get(eventduradj,'value');
eventdur=eventdur*10;
p=get(padj,'value');
lambda=get(lambdaadj,'value');
lambda=lambda*200000000;
limit0=get(startlimit,'value');
limit0=limit0*1000;
%def=get(default,'value');
if def==1
    lambda=100000000;
    p=0.6;
    thhold=5;
    eventdur=3;
    set(padj,'value',0.6);
    set(lambdaadj,'value',0.5);
     set(threshold,'value',0.5);
     set(eventduradj,'value',0.3);
      set(startlimit,'value',0.002);
end;
basel=baseline_(nt(num,:)',lambda,p);
%plot(basel,'r','LineWidth',1);
def=0;
% for i=0.1:0.01:0.29;
% plot(baseline_(nt(num,:)',lambda,p+i),'g','LineWidth',1);
% plot(baseline_(nt(num,:)',lambda,p-i),'g','LineWidth',1);
% end;
 j=1;
for i=1:length(imax)
    if basel(imax(i))<ymax(i)
        if abs(ymax(i)-basel(imax(i)))<10
noisemax(j)=abs(ymax(i)-basel(imax(i)));
j=j+1;
    end;
    end;
end;

j=1;
for i=1:length(imin)
    if basel(imin(i))>ymin(i)
        if abs(basel(imin(i))-ymin(i))<10
noisemin(j)=(abs(basel(imin(i))-ymin(i)));
j=j+1;
    end;
    end;
end;


noise=(mean(noisemin)+mean(noisemax))/2;
noisesd=(std(noisemin)+std(noisemax))/2;
% plot(basel+noise/2+noisesd,'r');
% plot(basel-noise/2-noisesd,'g');
% plot(basel+noise/2+noisesd+5,'b');
% plot(basel-noise/2-noisesd-5,'g');
clear('noisemin','noisemax');


smtlb = nt(num,:);%sgolayfilt(nt(num,:),21,41);   % 8 is a polynomial order.Apply the 3rd-order filter. out a significant portion of the signal's high frequency content along with the noise. Although Savitzky-Golay filters are more effective at preserving the pertinent high frequency components of the signal, they are less successful than standard averaging FIR filters at rejecting noise.
   %                                     %Savitzky-Golay filters are optimal in the sense that they minimize the least-squares error in fitting a polynomial to frames of noisy data.
 %plot(smtlb,'r','LineWidth',1);

%[ymax,imax,ymin,imin] = extrema(nt(num,:)); % num is the number of the actual cell cell N

% for i=1:length(imax)
%     if ymax(i)>basel(imax(i))+noise+noisesd+5
%      %   plot (imax(i),ymax(i),'r*')
%     end
% end
%      
% for i=1:length(imin)
%     if ymin(i)<basel(imin(i))-noise-noisesd-5
%       %  plot (imin(i),ymin(i),'g*')
%     end
% end

highthreshold=0.02+thhold; %1+thhold noise+noisesd+thhold;  
lowthreshold=0.01; %0.5 noise;  
plot(basel-highthreshold,'r');
plot(basel-lowthreshold,'g');

%Schmitt trigger algorithm for event detection
clear('indoff','inon','on','off');

onprox=1;
n=1;
i=1;
l=1;
while i<length(smtlb)
    i=i+n;
    n=1;
     
        if smtlb(i)<=basel(i)-highthreshold && i<length(nt(1,:)) % changed >+
        j=i;    %"event identification point"
        k=i;
        
            while smtlb(k)< basel(i)-lowthreshold && k<length(nt(1,:)) % changed >+
                  k=k+1;              
            end; 
       if k==length(nt(1,:))
           break;
       end;
             n=k-i;
             if j<=2 
                 j=3;
             end;
             if n>eventdur %4<n<30
            %onsets and offsets
             off(l)=min(nt(num,j-2:k));
             temp=find(nt(num,j-2:k)==min(nt(num,j-2:k))); %in case of more that one min found
             indoff(l)=temp(1);
             
             k=1;
             while imax(k)<j 
                 k=k+1;
             end;
           %  display(j);
           %  display(k-1);
          %   display(imax(k-1));
             %onprox=1;
             k=k-1;
            %proximity interval from "event identification point" (onprox) identification
            
             while max(nt(num,j-onprox:j))<basel(j-onprox) && j-onprox>1
                  onprox=onprox+1;
             end;  
           ii=k;
           %=========  
           
           %  plot (j,nt(num,j),'r*')
            k=indoff(l)+j-2-1 %angle algorithm
             if k<=limit0
             on(l)=nt(num,k);
              indon(l)=k;    
             end;    
            if k>limit0
            pa=abs(nt(num,k)-nt(num,k-1)); % (peak a) pa is a y catet
            pc=0.137; %pc is a x catet
            pb=sqrt(pa^2+pc^2); %peak hypotenuse
                
            angle0=acos((pa^2+pb^2-pc^2)/(2*pa*pb)); %angle between pa and pb
            angle=0; 
            
           
            % while ymax(k)<basel(imax(k))-lowthreshold && imax(k)>j-onprox && j-onprox>1 %
           % while nt(num,k)<basel(k)-lowthreshold && nt(num,k)>j-onprox && j-onprox>1 %
            while true 
            %  while angle<=angle0+0.05 && nt(num,k-1)>=nt(num,k) && j-onprox>1 %+0.05
                k=k-1;
                
                  pa=abs(nt(num,k)-nt(num,k-1)); %angle
                 pc=0.137; %angle
                 pb=sqrt(pa^2+pc^2);  %angle
                
                if nt(num,k-1)>=nt(num,k) 
                     angle=acos((pa^2+pb^2-pc^2)/(2*pa*pb)); %angle
                end;
                if nt(num,k-1)<nt(num,k) 
                  angle=3.14/2+acos((pa^2+pb^2-pc^2)/(2*pa*pb)); %angle
                end;
                
                
                 if angle<angle0
                   angle0=angle;  
                 end;
                 
                 pa0=abs(nt(num,k-1)-nt(num,k-2)) %angle
                 pc0=0.137; %angle
                 pb0=sqrt(pa0^2+pc0^2);  %angle
                if nt(num,k-2)>=nt(num,k-1)
                     angle00=acos((pa0^2+pb0^2-pc0^2)/(2*pa0*pb0));
                end;
                if nt(num,k-2)<nt(num,k-1)
                     angle00=3.14/2+acos((pa0^2+pb0^2-pc0^2)/(2*pa0*pb0));
                end; 
                
                 if angle00<=angle0 
                     k=k-1;
                     angle=angle00;
                 end;
          %   end;
          %   end; 
            if nt(num,k)>=basel(k)-lowthreshold-1 && angle>angle0+0.05
                break;
            end;
%             if nt(num,k)<basel(k-1)-lowthreshold
%                 k=k-1;
%             end;   
            end;
             on(l)=nt(num,k);
             indon(l)=k;
            % on(l)=ymax(k);
            % indon(l)=imax(k);
           end;
           if k==0 
               indoff(l)=1000;
               indon(l)=0;
           end;
             
%              on=max(nt(num,j-onprox:j));
%              indon=find(nt(num,j-onprox:j)==max(nt(num,j-onprox:j)));
           if  indoff(l)+j-1-indon(l)<400
               indoff(l)=indoff(l)+j-2-1;
                plot (j,smtlb(j),'b*')
                plot (indoff(l), off(l),'r.')
                offsetmat(num,indoff(l))=1;
                if l>1 && indon(l)==indon(l-1)
                    if ii>0
                     while ymax(ii)<=basel(imax(ii))-5 && imax(ii)>j-onprox && j-onprox>1
                       ii=ii-1;
                      end;
                       
                    
                    indon(l)=imax(ii);
                       on(l)=ymax(ii);
                                   
                    end;
                plot (indon(l), on(l), 'b.');
                onsetmat(num,indon(l-1))=1;
                offsetmat(num,indoff(l-1))=1;
                else
                 plot (indon(l), on(l), 'g.');
                 onsetmat(num,indon(l))=1;
                 offsetmat(num,indoff(l))=1;
               end;
                display(indoff(l));
                display(indon(l));
                l=l+1
           end;
             end; %if n>4
             if n>30
     
             end; %if n>30
        
        end; %if smtlb
       
end; %while
spk(num,:)=onsetmat(num,:);
dec(num,:)=offsetmat(num,:);
%end; %end of global for  
%===========================================================
%END OF A EVENT IDENTIFICATION ALGORITHM