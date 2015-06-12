% Frequency Response Analysis of Data
% 
% created by Srinivas Gorur-Shandilya at 3:19 , 10 June 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end
tic
warning off

%% Frequency Response Analysis
% In this document, we analyse data using frequency response plots (Bode plots, etc.)

%% Synthetic Data: Exponential filter with no noise
% In this section, we generate some synthetic data using some sine waves, a simple filter and no Gaussian additive noise:

t = 1e-3:1e-3:60;
a = randn(length(t),1) + .4*sin(t)';
K = filter_exp(30,1,1:500);
b = filter(K,1,a) + 0*randn(60e3,1);

figure('outerposition',[0 0 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on
subplot(2,8,1:6), hold on
ylabel('Input')
plot(t,a,'k')
subplot(2,8,9:14), hold on
plot(t,b,'r')
ylabel('Output')
subplot(2,8,7:8), hold on
plot(1e-3*(1:length(K)),K,'r')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end


%% 
% We now compute the frequency response by dividing the cross spectrum by the power spectrum, and compare it to the Fourier transform of the filter. We also use the NLID toolbox to calculate the same thing. 

% do your own compute
fs = 1e3;
[H,f] = tfestimate(a,b,[],[],[],fs);
ff = (0:length(K)-1)*(fs/length(K));
Kf = fft(K);

% use the toolbox
d = nldat;
set(d,'Data',[a b]);
set(d,'DomainValues',t);
set(d,'DomainIncr',1/fs);
set(d,'ChanNames',{'x','y'},'ChanUnits',{'x','y'})
FR = fresp(d);
fr = FR.data;
ft = 0:length(fr)-1;
ft = ft*get(FR,'DomainIncr');

figure('outerposition',[0 0 700 800],'PaperUnits','points','PaperSize',[700 800]); hold on
subplot(3,1,1), hold on
plot(f,abs(H))
plot(ft,abs(fr(:,1)))
plot(ff,abs(Kf))

legend({'abs(H)','toolbox result','FFT(K)'})
ylabel('Magnitude')
set(gca,'YScale','log','XLim',[0 500])

subplot(3,1,2), hold on
plot(f,-unwrap(rad2deg(angle(H))))
plot(ft,unwrap(rad2deg(angle(fr(:,1)))));
legend({'angle(H)','toolbox result'})
ylabel('Phase (degrees)')
set(gca,'YScale','linear','XLim',[0 500])

subplot(3,1,3), hold on
c=(abs(cpsd(a,b)).^2)./(cpsd(a,a).*cpsd(b,b));
plot(f,c)
plot(ft,fr(:,2))
legend({'Coherence','toolbox result'})
set(gca,'YLim',[0 1.1],'XLim',[0 500])
ylabel('Coherence')
xlabel('Frequency (Hz)')
PrettyFig;

if being_published
	snapnow
	delete(gcf)
end


%% 
% Can we get back the filter by inverse FTing the periodogram?

Khat = flipud(real(ifft(H)));
Khat = interp1(1:length(Khat),Khat,0.5:0.5:length(Khat));

Khat2 = revCorrFilter(a,b,'filter_length',500,'reg',0);

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(K,'k')
set(gca,'XLim',[0 500])
plot(Khat,'r')
plot(Khat2,'b')
xlabel('Filter lag (ms)')
legend('Actual Filter','IFFT(H)','Rev.Corr. Filter')
PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%%
% It's nice to see that both methods exactly get back the original filter. 

%% Synthetic Data: Exponential Filter with Output Noise
% In this section, we add some additive Gaussian Noise to the output and repeat the analysis. We shorten the nfft window to account for the noise, so we average over more Welch windows. 

t = 1e-3:1e-3:60;
a = randn(length(t),1) + .4*sin(t)';
K = filter_exp(30,1,1:500);
b = filter(K,1,a) + .1*randn(60e3,1);

fs = 1e3;
nfft = 2^10;
[H,f] = tfestimate(a,b,[],[],nfft,fs);
ff = (0:length(K)-1)*(1e3/length(K));
Kf = fft(K);

% use the toolbox
d = nldat;
set(d,'Data',[a b]);
set(d,'DomainValues',t);
set(d,'DomainIncr',1/fs);
set(d,'ChanNames',{'x','y'},'ChanUnits',{'x','y'})
FR = fresp(d);
fr = FR.data;
ft = 0:length(fr)-1;
ft = ft*get(FR,'DomainIncr');

figure('outerposition',[0 0 700 800],'PaperUnits','points','PaperSize',[700 800]); hold on
subplot(3,1,1), hold on
plot(ft,abs(fr(:,1)))
plot(f,abs(H))
plot(ff,abs(Kf))

legend({'toolbox result','abs(H)','FFT(K)'})
ylabel('Magnitude')
set(gca,'YScale','log','XLim',[0 500],'YLim',[.01 2])

subplot(3,1,2), hold on
plot(f,unwrap(rad2deg(angle(H))))
plot(ft,unwrap(rad2deg(angle(fr(:,1)))));
legend({'angle(H)','toolbox result'})
ylabel('Phase (degrees)')
set(gca,'YScale','linear','XLim',[0 500])

subplot(3,1,3), hold on
[c,f] = mscohere(a,b,[],[],[],fs);
plot(ft,fr(:,2))
plot(f(1:10:end),csaps(f,c(:,end),1e-3,f(1:10:end)));
legend({'toolbox result','Smoothed Coherence'})
set(gca,'YLim',[0 1.1],'XLim',[0 500])
ylabel('Coherence')
xlabel('Frequency (Hz)')
PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%% 
% This looks pretty bad. Can we extract filters? How bad does the FFT do?


Khat = (real(ifft(H)));
Khat = interp1(1:length(Khat),Khat,0.5:0.5:length(Khat));

K_toolbox = (real(ifft(fr(:,1))));
K_toolbox = interp1(1:length(K_toolbox),K_toolbox,0.5:0.5:length(K_toolbox));

Khat2 = revCorrFilter(a,b,'filter_length',500,'reg',.01);

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(K,'k')
set(gca,'XLim',[0 500])
plot(Khat,'r')
plot(Khat2,'b')
plot(K_toolbox,'g')
xlabel('Filter lag (ms)')
legend('Actual Filter','IFFT(H)','Rev.Corr. Filter','IFFT(toolbox)')
PrettyFig;

if being_published
	snapnow
	delete(gcf)
end


%%
% So even though the filter extraction is really good (indicating the noise is weak), the Bode diagrams are noisy, and it's hard to see anything there. 


%% Synthetic Data: Ankle Joint Compliance Data
% This data is synthesised by some models used in Westwick and Kearney. We're going to attempt recreate Fig 5.7 (pg 118). 

iodata = ExampleData(1,10); d = iodata.Data; a = d(:,1); b = d(:,2);
t = 0.002*(1:length(a));

fs = 1/0.002;
nfft = 2^8;
[H,f] = tfestimate(a,b,[],[],nfft,fs);

% use the toolbox
FR = fresp(iodata);
fr = FR.data;
ft = 0:length(fr)-1;
ft = ft*get(FR,'DomainIncr');

figure('outerposition',[0 0 700 800],'PaperUnits','points','PaperSize',[700 800]); hold on
subplot(3,1,1), hold on
plot(f,abs(H))
plot(ft,abs(fr(:,1)))
legend({'abs(H)','toolbox result'})
ylabel('Magnitude')
set(gca,'YScale','log','XLim',[0 250])

subplot(3,1,2), hold on
plot(f,unwrap(rad2deg(angle(H))))
plot(ft,unwrap(rad2deg(angle(fr(:,1)))));
legend({'angle(H)','toolbox result'})
ylabel('Phase (degrees)')
set(gca,'YScale','linear','XLim',[0 250])

subplot(3,1,3), hold on
[c,f] = mscohere(a,b,[],[],[],fs);
plot(ft,fr(:,2))
plot(f(1:10:end),csaps(f,c(:,end),1e-3,f(1:10:end)));
legend({'toolbox result','Smoothed coherence'})
set(gca,'YLim',[0 1.1],'XLim',[0 250])
ylabel('Coherence')
xlabel('Frequency (Hz)')
PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%% Real Data
% We now analyse real data: flickering stimulus presented to a couple of ab3A neurons 20 times. We plot the real and imaginary parts of the transfer functions, with the coherence on a trial-wise basis.  


load('/local-data/DA-paper/large-variance-flicker/2015_01_28_CS_ab3_2_EA.mat')
PID = data(4).PID;
time = 1e-4*(1:length(PID));
all_spikes = spikes(4).A;
B_spikes = spikes(4).B;
load('/local-data/DA-paper/large-variance-flicker/2015_01_28_CS_ab3_3_EA.mat')
PID = vertcat(PID,data(4).PID);
all_spikes = vertcat(all_spikes,spikes(4).A);
B_spikes = vertcat(B_spikes,spikes(4).B);

% A spikes --> firing rate
hash = DataHash(full(all_spikes));
cached_data = cache(hash);
if isempty(cached_data)
	fA = spiketimes2f(all_spikes,time);
	cache(hash,fA);
else
	fA = cached_data;
end

tA = 1e-3*(1:length(fA));
PID2 = fA;
for i = 1:width(PID2)
	PID2(:,i) = interp1(time,PID(i,:),tA);
end
PID = PID2; clear PID2
% some minor cleaning up
PID(end,:) = PID(end-1,:); 


% compute transfer function on a trial-wise basis
fs = 1e3;
clear H
for i = 1:20
	a = PID(:,i);
	b  =fA(:,i);
	a = a(20e3:end);
	b = b(20e3:end);
	t = 1e-3*(1:length(a));
	[H(i,:),f] = tfestimate(a,b,[],[],2^10,fs);

	% also use the toolbox
	d = nldat;
	set(d,'Data',[a b]);
	set(d,'DomainValues',t);
	set(d,'DomainIncr',1/fs);
	set(d,'ChanNames',{'x','y'},'ChanUnits',{'x','y'})
	FR = fresp(d);
	fr = FR.data;
	ft = 0:length(fr)-1;
	ft = ft*get(FR,'DomainIncr');
	% save it
	G(i,:) = fr(:,1);
	C(i,:) = fr(:,2);

end

figure('outerposition',[0 0 900 700],'PaperUnits','points','PaperSize',[900 700]); hold on
subplot(3,1,1), hold on
A = f*0;
for i = 1:20
	A = A + abs(H(i,:))';
end
plot(f,A/20,'k')
set(gca,'XScale','log')

A = ft*0;
for i = 1:20
	A = A + abs(G(i,:));
end
plot(ft,A/20,'r')
legend('Abs(H)','toolbox')
ylabel('Amplitude')

subplot(3,1,2), hold on
A = f*0;
for i = 1:20
	A = A + angle(H(i,:))';
end
plot(f,rad2deg(A/20),'k')
A = ft*0;
for i = 1:20
	A = A + angle(G(i,:));
end
plot(ft,rad2deg(A/20),'r')
legend('Angle(H)','toolbox')
set(gca,'XScale','log')
ylabel('Amplitude')

set(gca,'XScale','log','YLim',[-180 180])
ylabel('Phase (degrees)')

c = [];
subplot(3,1,3), hold on
for i = 1:20
	a = PID(:,i);
	b  =fA(:,i);
	a = a(20e3:end);
	b = b(20e3:end);

	[temp,f] = mscohere(a,b,[],[],[],fs);
	c(i,:) = csaps(f,temp,.5,f(1:10:end));
end
plot(f(1:10:end),mean2(c),'k')
plot(ft,mean2(C),'r')
set(gca,'XScale','log')
legend('Smoothed Coherence','toolbox')
ylabel('Coherence')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end


%% Version Info
% The file that generated this document is called:
disp(mfilename)

%%
% and its md5 hash is:
Opt.Input = 'file';
disp(DataHash(strcat(mfilename,'.m'),Opt))

%%
% This file should be in this commit:
[status,m]=unix('git rev-parse HEAD');
if ~status
	disp(m)
end

t = toc;
%% 
% This document was built in: 
disp(strcat(oval(t,3),' seconds.'))


