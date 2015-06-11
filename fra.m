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
% We now compute the frequency response by dividing the cross spectrum by the power spectrum, and compare it to the Fourier transform of the filter:

fs = 1e3;
[~,f] = pwelch(a,[],[],[],fs);
H = cpsd(a,b)./cpsd(a,a);
ff = (0:length(K)-1)*(1e3/length(K));
Kf = fft(K);

figure('outerposition',[0 0 700 800],'PaperUnits','points','PaperSize',[700 800]); hold on
subplot(3,1,1), hold on
plot(ff,abs(Kf))
plot(f,abs(H))
legend({'real(H)','FFT(K)'})
ylabel('Magnitude')
set(gca,'YScale','log','XLim',[0 500])

subplot(3,1,2), hold on
plot(f,unwrap(rad2deg(angle(H))))
ylabel('Phase (degrees)')
set(gca,'YScale','linear','XLim',[0 500])


subplot(3,1,3), hold on
c=(abs(cpsd(a,b)).^2)./(cpsd(a,a).*cpsd(b,b));
plot(f,c)
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
% In this section, we add some additive Gaussian Noise to the output and repeat the analysis. 
t = 1e-3:1e-3:60;
a = randn(length(t),1) + .4*sin(t)';
K = filter_exp(30,1,1:500);
b = filter(K,1,a) + .1*randn(60e3,1);

fs = 1e3;
[~,f] = pwelch(a,[],[],[],fs);
H = cpsd(a,b)./cpsd(a,a);
ff = (0:length(K)-1)*(1e3/length(K));
Kf = fft(K);

figure('outerposition',[0 0 700 800],'PaperUnits','points','PaperSize',[700 800]); hold on
subplot(3,1,1), hold on
plot(ff,abs(Kf))
plot(f,abs(H))
legend({'real(H)','FFT(K)'})
ylabel('Magnitude')
set(gca,'YScale','log','XLim',[0 500])

subplot(3,1,2), hold on
plot(f,unwrap(rad2deg(angle(H))))
ylabel('Phase (degrees)')
set(gca,'YScale','linear','XLim',[0 500])


subplot(3,1,3), hold on
c=(abs(cpsd(a,b)).^2)./(cpsd(a,a).*cpsd(b,b));
plot(f,c)
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


Khat = flipud(real(ifft(H)));
Khat = interp1(1:length(Khat),Khat,0.5:0.5:length(Khat));

Khat2 = revCorrFilter(a,b,'filter_length',500,'reg',.01);

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
% So even though the filter extraction is really good (indicating the noise is weak), the Bode diagrams are rubbish, and it's hard to see anything there. 


%% Synthetic Data: Ankle Joint Compliance Data
% This data is synthesised by some models used in Westwick and Kearney. We're going to attempt recreate Fig 5.7 (pg 118). 

iodata = ExampleData(1,10); d = iodata.Data; a = d(:,1); b = d(:,2);
t = 0.002*(1:length(a));

fs = 1/0.002;
[denom,f] = pwelch(a,[],[],[],fs);
H = cpsd(a,b)./cpsd(a,a);

figure('outerposition',[0 0 700 800],'PaperUnits','points','PaperSize',[700 800]); hold on
subplot(3,1,1), hold on
plot(f,abs(H))
ylabel('Magnitude')
set(gca,'YScale','log')

subplot(3,1,2), hold on
plot(f,unwrap(rad2deg(angle(H))))
ylabel('Phase (degrees)')

subplot(3,1,3), hold on
c=(abs(cpsd(a,b)).^2)./(cpsd(a,a).*cpsd(b,b));
plot(f,c)
ylabel('Coherence')
PrettyFig;

if being_published
	snapnow
	delete(gcf)
end



return

%% Phase Relationships and Coherence
% In the previous section, we plotted the real part of the transfer function H. Now, we plot the imaginary part, which tells us how the phase depends on the frequency. 

figure('outerposition',[0 0 800 800],'PaperUnits','points','PaperSize',[800 800]); hold on
subplot(2,1,1), hold on
plot(f,imag(H));
ylabel('Phase (degrees)')
xlabel('Frequency (Hz)')
set(gca,'XScale','log')

subplot(2,1,2), hold on
c=(abs(cpsd(a,b)).^2)./(pwelch(a).*pwelch(b));
plot(f,c)
ylabel('Coherence')
xlabel('Frequency (Hz)')
set(gca,'XScale','log')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%% More synthetic data: white noise analysis
% Now we repeat the analysis, but with white noise. 


t = 1e-3:1e-3:60;
a = randn(length(t),1); 
K = filter_alpha2(30,70,1,.2,1:500);
b = filter(K,1,a) + randn(60e3,1);

[denom,f] = pwelch(a,[],[],[],1e3);
H = cpsd(a,b)./denom;
ff = (0:499)*(1e3/500);
Kf = fft(K);

figure('outerposition',[0 0 700 900],'PaperUnits','points','PaperSize',[700 900]); hold on
subplot(3,1,1), hold on
plot(f,real(H)/max(real(H)),'r')
plot(ff,real(Kf)/max(real(Kf)));
legend({'real(H)','FFT(K)'})
ylabel('Power (norm)')
xlabel('Frequency (Hz)')
set(gca,'XScale','log')

subplot(3,1,2), hold on
plot(f,imag(H));
ylabel('Phase (degrees)')
xlabel('Frequency (Hz)')
set(gca,'XScale','log')

subplot(3,1,3), hold on
c=(abs(cpsd(a,b)).^2)./(pwelch(a).*pwelch(b));
plot(f,c)
ylabel('Coherence')
xlabel('Frequency (Hz)')
set(gca,'XScale','log')

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
clear H
for i = 1:20
	a = PID(:,i);
	b  =fA(:,i);
	a = a(20e3:end);
	b = b(20e3:end);
	[denom,f] = pwelch(a,[],[],[],1e3);
	H(i,:) = cpsd(a,b)./denom;
end

figure('outerposition',[0 0 500 900],'PaperUnits','points','PaperSize',[500 900]); hold on
subplot(3,1,1), hold on
cmap = parula(21);
for i = 1:20
	plot(f,reala(H(i,:)),'Color',cmap(i,:))
end
set(gca,'XScale','log')
ylabel('Amplitude')

subplot(3,1,2), hold on
cmap = parula(21);
for i = 1:20
	plot(f,imag(H(i,:)),'Color',cmap(i,:))
end
set(gca,'XScale','log','YLim',[-180 180],'XLim',[1e-3 1e3])
ylabel('Phase (degrees')

subplot(3,1,3), hold on
cmap = parula(21);
for i = 1:20
	a = PID(:,i);
	b  =fA(:,i);
	a = a(20e3:end);
	b = b(20e3:end);
	c=(abs(cpsd(a,b)).^2)./(pwelch(a).*pwelch(b));
	plot(f,c,'Color',cmap(i,:))
end
set(gca,'XScale','log')
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


