% Chrimson_MSG.m
% 
% created by Srinivas Gorur-Shandilya at 4:24 , 12 November 2015. Contact me at http://srinivas.gs/contact/
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

%% Gain changes with light in Chrimson flies
% In this document we look at how gain changes using changes in light mean and variance in ORNs expressing Chrimson. 

%% Gain changes with mean light. 
% In this section we look at how ORN gain changes with light background, where we measure the gain in response to a fluctuating light stimulus. No odour at all in this experiment. 

p = '/local-data/DA-paper/Chrimson/MSG/ok';
[~, ~, fA, paradigm, orn, fly, AllControlParadigms] = consolidateData(p,true);


% make new vectors for mean and range of the stimulus
LED = NaN*fA;
r = NaN*paradigm; m = NaN*paradigm;
for i = 1:length(paradigm)
	temp = AllControlParadigms(paradigm(i)).Name;
	r(i) = str2double(temp(3:strfind(temp,'_')-1));
	m(i) = str2double(temp(strfind(temp,'_')+3:end));
	LED(:,i) = AllControlParadigms(paradigm(i)).Outputs(1,1:10:end);
end

% extract filters
a = 25e3; z = 55e3;
[K,fp,gain,gain_err] = extractFilters(LED,fA,'use_cache',true,'a',a,'z',z);

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
nplots = length(unique(r));
all_r = unique(r);
for i = 1:nplots
	autoPlot(nplots,i); hold on
	plot_this = r == all_r(i);
	mean_stim = m(plot_this);
	g = gain(plot_this);
	plot(mean_stim,g,'k+')
	title(['Range=' oval(all_r(i))])
	xlabel('Mean Light (V)')
	ylabel('Gain (Hz/V)')
end

prettyFig()

if being_published
	snapnow
	delete(gcf)
end

%% Gain Changes with Light contrast 
% Now we look at how gain changes with the contrast of the light. 

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(r(m==2),gain(m==2),'k+')
xlabel('Range of Light Input (V)')
ylabel('Gain (Hz/V)')
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
[status,git_hash]=unix('git rev-parse HEAD');
if ~status
	disp(git_hash)
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
