% STC_Analayis.m
%
% created by Srinivas Gorur-Shandilya at 11:08 , 02 November 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.



% add homebrew path
path1 = getenv('PATH');
if isempty(strfind(path1,':/usr/local/bin'))
    path1 = [path1 ':/usr/local/bin'];
end
setenv('PATH', path1);

% this code determines if this function is being called
calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
		unix(['tag -a publish-failed ',which(mfilename)]);
		unix(['tag -r published ',which(mfilename)]);
	end
end
tic


%% STC Analysis
% In this document we analyse the spike-triggered covariance using some synthetic data. Here, we generate stimuli using some random Gaussian noise, use an exponential filter and a threshold to generate spikes, and then look at the STC of that. 

RandStream.setGlobalStream(RandStream('mt19937ar','Seed',1984)); 
s = randn(60e3,1);
K = filter_exp(10,1,1:1000);
K = K/max(K);
r = filter(K,1,s);
r(r < mean(r)+std(r)) = 0;
r(r>0) = 1;

[eigen_vectors, eigen_values,C_prior,C_spike] = STC(r,s);

figure('outerposition',[0 0 1200 900],'PaperUnits','points','PaperSize',[1200 900]); hold on
subplot(3,4,1:3), hold on
plot(1e-4*(1:length(s)),s)
set(gca,'XLim',[1 1.1])
ylabel('Stimulus')

subplot(3,4,4), hold on
imagesc(C_prior)
title('C_{prior}')
axis tight
colorbar

subplot(3,4,5:7), hold on
raster2(r)
title('Spikes')
set(gca,'XLim',[1 1.1])

subplot(3,4,8), hold on
imagesc(C_spike)
title('C_{spike}')
axis tight
colorbar

subplot(3,2,5), hold on
l(1) = plot(K,'k');
Khat = eigen_vectors(:,1);
Khat = Khat - mean(Khat(end-100:end));
Khat = Khat/max(Khat);
l(2) = plot(Khat,'r');
legend(l,{'Actual Filter','Rescaled first Eigenvector'})
set(gca,'XLim',[0 500])

subplot(3,2,6), hold on
plot(eigen_values,'k-+')
title('First six Eigenvalues')

prettyFig()

if being_published
	snapnow
	delete(gcf)
end

%% A more complicated filter
% Now we repeat the analysis with a more realistic filter:

s = randn(60e3,1);
p.tau1 = 10;
p.tau2 = 30;
p.A = 2;
p.n = 2;
K = filter_gamma2(1:500,p);
K = K/max(K);
r = filter(K,1,s);
r(r < mean(r)+std(r)) = 0;
r(r>0) = 1;

[eigen_vectors, eigen_values,C_prior,C_spike] = STC(r,s);

figure('outerposition',[0 0 1200 900],'PaperUnits','points','PaperSize',[1200 900]); hold on
subplot(3,4,1:3), hold on
plot(1e-4*(1:length(s)),s)
set(gca,'XLim',[1 1.1])
ylabel('Stimulus')

subplot(3,4,4), hold on
imagesc(C_prior)
title('C_{prior}')
axis tight
colorbar

subplot(3,4,5:7), hold on
raster2(r)
title('Spikes')
set(gca,'XLim',[1 1.1])

subplot(3,4,8), hold on
imagesc(C_spike)
title('C_{spike}')
axis tight
colorbar

subplot(3,2,5), hold on
l(1) = plot(K,'k');
Khat = eigen_vectors(:,1);
Khat = Khat - mean(Khat(end-100:end));
Khat = Khat/max(Khat);
l(2) = plot(Khat,'r');
legend(l,{'Actual Filter','Rescaled first Eigenvector'})
set(gca,'XLim',[0 500])

subplot(3,2,6), hold on
plot(eigen_values,'k-+')
title('First six Eigenvalues')

prettyFig()

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
disp(dataHash(strcat(mfilename,'.m'),Opt))

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

% tag the file as being published 

if being_published
	unix(['tag -a published ',which(mfilename)]);
	unix(['tag -r publish-failed ',which(mfilename)]);
end
