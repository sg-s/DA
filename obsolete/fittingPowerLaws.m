% fittingPowerLaws.m
% this document looks at how we can fit power laws to data, and the errors that arise from various ways of doing things.
% 
% created by Srinivas Gorur-Shandilya at 2:57 , 06 July 2015. Contact me at http://srinivas.gs/contact/
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

%% How do absolute errors in measurement affect our ability to fit a power law to the data?
% Imagine that we obtain a few samples from a function that varies as a inverse power low. However, the samples have additive Gaussian noise on them. How will our ability to fit a line vary with the amount of Gaussian noise we add to it?

a = 100;
f = @(x) (a./x);
noise = linspace(0,5,3);
n = 50;

figure('outerposition',[0 0 1300 500],'PaperUnits','points','PaperSize',[1300 500]); hold on 
for i = 1:length(noise)
	subplot(1,length(noise),i), hold on
	x = 1e-3:1000;
	plot(x,f(x),'k--')
	x = rand(1,100).*logspace(-1,2,100);
	y = f(x) + noise(i)*randn(1,100);
	plot(x,y,'k+')

	% now fit a power law
	rm_this = x<0 | y <0;
	y(rm_this) = [];
	x(rm_this) = [];
	ff = fit(x(:),y(:),'power1');

	plot(x,ff(x),'r')
	xlabel('x')
	ylabel('f(x)')

	title(strcat('b=',oval(ff.b)))
	set(gca,'XScale','log','YScale','log','XLim',[.1 100],'YLim',[1e-1 1e3])
end


PrettyFig

if being_published
	snapnow
	delete(gcf)
end

%%
% Now we repeat this hundreds of times, to get a estimate of the spread

noise = linspace(0,5,30);
nrep = 30;
b_mean = [];
b_std = [];

for i = 1:length(noise)
	textbar(i,length(noise))
	temp = [];
	for j = 1:nrep
		x = rand(1,100).*logspace(-1,2,100);
		y = f(x) + noise(i)*randn(1,100);

		% now fit a power law
		rm_this = x<0 | y <0;
		y(rm_this) = [];
		x(rm_this) = [];
		ff = fit(x(:),y(:),'power1');

		temp = [temp ff.b];
	end
	b_mean(i) = mean(temp);
	b_std(i) = std(temp);
end

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
errorbar(noise,(b_mean+1),b_std)
xlabel('Noise Level')
ylabel('Actual Slope - Estimated Slope')

PrettyFig

if being_published
	snapnow
	delete(gcf)
end

%%
% This shows that even large amounts of additive noise can't change the estimate of the exponent by much. 

%% Relative Errors
% Now, we repeat the simulation, but assume that the errors in the measurement are relative. 

noise = linspace(0,.5,30);
nrep = 30;
b_mean = [];
b_std = [];

for i = 1:length(noise)
	textbar(i,length(noise))
	temp = [];
	for j = 1:nrep
		x = rand(1,100).*logspace(-1,2,100);
		y = f(x);
		y = y + y.*noise(i)*randn(1,100)';

		% now fit a poywer law
		rm_this = x<0 | y <0;
		if sum(rm_this) < 90
			
			y(rm_this) = [];
			x(rm_this) = [];
			ff = fit(x(:),y(:),'power1');

			temp = [temp ff.b];
		end
	end
	b_mean(i) = mean(temp);
	b_std(i) = std(temp)/sqrt(length(temp));
end

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
errorbar(noise,(b_mean+1),b_std)
xlabel('Noise Level')
ylabel('Actual Slope - Estimated Slope')

PrettyFig

if being_published
	snapnow
	delete(gcf)
end

%%
% So there is a systematic error in estimating slope when errors are relative (which is natural). 

%% Fitting to the log transform
% What if we first take the log of the data and then fit a line to this? 

noise = linspace(0,.5,30);
nrep = 30;
b_mean = [];
b_std = [];

for i = 1:length(noise)
	textbar(i,length(noise))
	temp = [];
	for j = 1:nrep
		x = rand(1,100).*logspace(-1,2,100);
		y = f(x);
		y = y + y.*noise(i)*randn(1,100)';

		% now fit a power law
		rm_this = x<0 | y <0;
		if sum(rm_this) < 90
			
			y(rm_this) = [];
			x(rm_this) = [];
			ff = fit(log(x(:)),log(y(:)),'poly1');

			temp = [temp ff.p1];
		end
	end
	b_mean(i) = mean(temp);
	b_std(i) = std(temp)/sqrt(length(temp));
end

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
errorbar(noise,(b_mean+1),b_std)
xlabel('Noise Level')
ylabel('Actual Slope - Estimated Slope')

PrettyFig

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
