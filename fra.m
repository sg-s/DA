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

%% Synthetic Data
% In this section, we generate some synthetic data using some sine waves, a simple filter and some Gaussian additive noise:

t = 1e-3:1e-3:60;
a = randn(length(t),1) + .4*sin(t)';
K = filter_alpha2(30,70,1,.2,1:500);
b = filter(K,1,a) + randn(60e3,1);

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

[denom,f] = pwelch(a,[],[],[],1e3);
H = cpsd(a,b)./denom;
ff = (0:499)*(1e3/500);
Kf = fft(K);

figure('outerposition',[0 0 700 500],'PaperUnits','points','PaperSize',[700 500]); hold on
plot(f,real(H)/max(real(H)),'r')
plot(ff,real(Kf)/max(real(Kf)));
legend({'real(H)','FFT(K)'})
ylabel('Power (norm)')
xlabel('Frequency (Hz)')
set(gca,'XScale','log')
PrettyFig;

if being_published
	snapnow
	delete(gcf)
end


%% 
% Can we get back the filter by inverse FTing the periodogram?

Khat = flipud(real(ifft(H)));
Khat = Khat/max(Khat);
Khat = interp1(1:length(Khat),Khat,0.5:0.5:length(Khat));

[Khat2,~,filtertime] = FindBestFilter(a,b,[],'filter_length=600;','regmax = .01;','regmin=.01;');
Khat2 = Khat2/max(Khat2);

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(K,'k')
set(gca,'XLim',[0 500])
plot(Khat,'r')
plot(filtertime,Khat2,'b')
xlabel('Filter lag (ms)')
legend('Actual Filter','IFFT(H)','Rev.Corr. Filter')
PrettyFig;



if being_published
	snapnow
	delete(gcf)
end

%% Phase Relationships
% In the previous section, we plotted the real part of the transfer function H. Now, we plot the imaginary part, which tells us how the phase depends on the frequency. 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(f,imag(H));
ylabel('Phase (degrees)')
xlabel('Frequency (Hz)')
set(gca,'XScale','log')

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


