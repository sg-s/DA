function [spkvec,timevec]=spikegauss(timestamps,srate,min_timevec,max_timevec,alpha)
% Generate SPKVEC time series from TIMESTAMPS. 
% Each spike is represented by a gaussian with maximum value 1 on the
% timestamp. ALPHA is inversely proportional to standard deviation.
%
% If srate=1000, 95% of the gaussian is included in timestamp+-200ms when alpha=5, 
% 100ms when alpha=10, 50ms when alpha=20, 25ms when alpha=40.
%
% [spkvec,timevec]=spikegauss(timestamps,srate,min_timevec,max_timevec,alpha)
%
% Example :
%
%     timestamps=[-1.22 0.33 0.34 0.35 0.40 3.70 7.30]; %sec
%     srate=1000;
%     min_timevec=-4; 
%     max_timevec=8;
%     alpha=20;
%
%  [spkvec,timevec]=spikegauss(timestamps,srate,min_timevec,max_timevec,alpha);
%
%     plot(timevec,spkvec,'k')
%     hold on
%     plot([min_timevec max_timevec],[1 1],'-r')
%     plot(timestamps,randn(size(timestamps))/10+1,'ob')
%     hold off
%
% Doubts, bugs: rpavao@gmail.com

% %%
% clc
% clear
% timestamps=[-1.22 0.33 0.34 0.35 0.40 3.70 7.30]; %sec
% srate=1000; %only even numbers...
% min_timevec=-4;
% max_timevec=8;
% alpha=20;

spkpos=round( (timestamps-min_timevec) * srate); %sec
timevec=0:1/srate:max_timevec-min_timevec;
spkvec=zeros(size(timevec));
    
for i=1:length(spkpos)
    temp=gausswin(srate-1,alpha)';
    start_end=[spkpos(i)-(round(srate/2)-1) spkpos(i)+(round(srate/2)-1) 1 srate-1];
    if spkpos(i)<=(round(srate/2)-1); start_end(1)=1; start_end(3)=(round(srate/2)+1)-spkpos(i); end
    if spkpos(i)+(round(srate/2)-1)>=length(spkvec); start_end(2)=length(spkvec); start_end(4)=1+start_end(2)-start_end(1); end    
    spkvec(start_end(1):start_end(2))=spkvec(start_end(1):start_end(2))+temp(start_end(3):start_end(4));
end

timevec=timevec+min_timevec+1/srate;

% plot(timevec,spkvec,'k')
% hold on
% plot([min_timevec max_timevec],[1 1],'-r')
% plot(timestamps,randn(size(timestamps))/10+1,'ob')
% hold off   

