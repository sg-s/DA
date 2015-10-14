% makeFig7.m
% makes fig 7 for the paper, showing gain in the LFP
% 
% created by Srinivas Gorur-Shandilya at 2:01 , 12 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


% add homebrew path
path1 = getenv('PATH');
if isempty(strfind(path1,[':/usr/local/bin']))
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

%       ##      ## ######## ########  ######## ########  
%       ##  ##  ## ##       ##     ## ##       ##     ## 
%       ##  ##  ## ##       ##     ## ##       ##     ## 
%       ##  ##  ## ######   ########  ######   ########  
%       ##  ##  ## ##       ##     ## ##       ##   ##   
%       ##  ##  ## ##       ##     ## ##       ##    ##  
%        ###  ###  ######## ########  ######## ##     ## 

%       ######## ########  ######  ##     ## ##    ## ######## ########  
%       ##       ##       ##    ## ##     ## ###   ## ##       ##     ## 
%       ##       ##       ##       ##     ## ####  ## ##       ##     ## 
%       ######   ######   ##       ######### ## ## ## ######   ########  
%       ##       ##       ##       ##     ## ##  #### ##       ##   ##   
%       ##       ##       ##    ## ##     ## ##   ### ##       ##    ##  
%       ##       ########  ######  ##     ## ##    ## ######## ##     ## 


[PID, LFP, fA, paradigm, orn, fly, AllControlParadigms, paradigm_hashes] = consolidateData('/local-data/DA-paper/LFP-MSG/september',1);


% remove baseline from all PIDs
for i = 1:width(PID)
	PID(:,i) = PID(:,i) - mean(PID(1:5e3,i));
end

% sort the paradigms sensibly
sort_value = [];
for i = 1:length(AllControlParadigms)
	sort_value(i) = (mean(AllControlParadigms(i).Outputs(1,:)));
end
[~,idx] = sort(sort_value);
AllControlParadigms = AllControlParadigms(idx);
paradigm_new = paradigm*NaN;
for i = 1:length(idx)
	paradigm_new(paradigm == idx(i)) = i;
end
paradigm = paradigm_new;

% remove "Flicker" from paradigm names
for i = 1:length(AllControlParadigms)
	AllControlParadigms(i).Name = strrep(AllControlParadigms(i).Name,'Flicker-','');
end


% throw out trials where we didn't record the LFP, for whatever reason
not_LFP = find((max(abs(LFP))) < 0.1);
LFP(:,not_LFP) = NaN;

% throw our bad traces
bad_trials = (sum(fA) == 0 | isnan(sum(fA)) |  isnan(sum(LFP)));
LFP(:,bad_trials) = [];
PID(:,bad_trials) = [];
fA(:,bad_trials) = [];
paradigm(bad_trials) = [];
orn(bad_trials) = [];

% band pass all the LFP
filtered_LFP = LFP;
for i = 1:width(LFP)
	filtered_LFP(:,i) = 10*bandPass(LFP(:,i),1000,10); % now in mV
end

figure('outerposition',[0 0 1600 800],'PaperUnits','points','PaperSize',[1600 800]); hold on
axes_handles(1) = subplot(2,4,1:2);
axes_handles(2) = subplot(2,4,5);
axes_handles(3) = subplot(2,4,6);

axes_handles(4) = subplot(2,4,3);
axes_handles(5) = subplot(2,4,7);
axes_handles(6) = subplot(2,4,4);
axes_handles(7) = subplot(2,4,8);

for i = 1:length(axes_handles)
	hold(axes_handles(i),'on')
end

% show the filtered LFP
example_orn = 4;
these_paradigms = unique(paradigm(:,orn == example_orn));
c = parula(1+length(these_paradigms));
time = 1e-3*(1:length(LFP(40e3:55e3,1))) + 40;
for i = 1:length(these_paradigms)
	plot_this = find(orn == example_orn & paradigm == these_paradigms(i));
	for j = 1:length(plot_this)
		plot(axes_handles(1),time,filtered_LFP(40e3:55e3,plot_this(j)),'Color',c(i,:))
	end
end

xlabel(axes_handles(1),'Time (s)')
ylabel(axes_handles(1),'\DeltaLFP (mV)')

%     ##       ######## ########      ######      ###    #### ##    ## 
%     ##       ##       ##     ##    ##    ##    ## ##    ##  ###   ## 
%     ##       ##       ##     ##    ##         ##   ##   ##  ####  ## 
%     ##       ######   ########     ##   #### ##     ##  ##  ## ## ## 
%     ##       ##       ##           ##    ##  #########  ##  ##  #### 
%     ##       ##       ##           ##    ##  ##     ##  ##  ##   ### 
%     ######## ##       ##            ######   ##     ## #### ##    ## 


% show the lfp gain with mean stimulus
a = 10e3; z = 50e3;
[K1,LFP_pred,LFP_gain,LFP_gain_err] = extractFilters(PID,filtered_LFP,'use_cache',true,'a',a,'z',z);

ss  =50;
c= parula(max(paradigm)+1);
l = [];
filtertime = 1e-3*(1:length(K1))-.1;

for i = 1:width(LFP)
	plot(axes_handles(2),mean(PID(a:z,i)),LFP_gain(i),'+','Color',c(paradigm(i),:))
end


x = mean(PID(a:z,:)); x = x(:);
y = LFP_gain(:);

xx = NaN(length(unique(paradigm)),1);
yy = NaN(length(unique(paradigm)),1);
ww = NaN(length(unique(paradigm)),1);

for i = 1:length(yy)
	xx(i) = mean(x(paradigm==i));
	yy(i) = mean(y(paradigm==i));
	ww(i) = 1./sem(y(paradigm==i));
end


ff = fit(xx(:),yy(:),'power1','Weights',ww);
clear l

l(1) = plot(axes_handles(2),sort(mean(PID(a:z,:))),ff(sort(mean(PID(a:z,:)))),'k--');
L = {};
L{1} = ['y = \alpha x^{', oval(ff.b), '}'];

fo = fitoptions('power1');
fo.Upper = [NaN -1];
fo.Lower = [NaN -1];
fo.Weights = ww;
ff = fit(xx(:),yy(:),'power1',fo);
l(2) = plot(axes_handles(2),sort(mean(PID(a:z,:))),ff(sort(mean(PID(a:z,:)))),'k');
L{2} = ['y = \alpha x^{', oval(ff.b), '}'];

legend(l,L)
legend('boxoff')
set(axes_handles(2),'XScale','log','YScale','log')
xlabel(axes_handles(2),'Mean Stimulus (V)')
ylabel(axes_handles(2),'LFP Gain (mV/V)')

% ######## #### ########  #### ##    ##  ######       ######      ###    #### ##    ## 
% ##        ##  ##     ##  ##  ###   ## ##    ##     ##    ##    ## ##    ##  ###   ## 
% ##        ##  ##     ##  ##  ####  ## ##           ##         ##   ##   ##  ####  ## 
% ######    ##  ########   ##  ## ## ## ##   ####    ##   #### ##     ##  ##  ## ## ## 
% ##        ##  ##   ##    ##  ##  #### ##    ##     ##    ##  #########  ##  ##  #### 
% ##        ##  ##    ##   ##  ##   ### ##    ##     ##    ##  ##     ##  ##  ##   ### 
% ##       #### ##     ## #### ##    ##  ######       ######   ##     ## #### ##    ## 


a = 10e3; z = 50e3;
[K2,fA_pred,fA_gain,fA_gain_err] = extractFilters(PID,fA,'use_cache',true,'a',a,'z',z);

for i = 1:width(fA)
	plot(axes_handles(3),mean(PID(a:z,i)),fA_gain(i),'+','Color',c(paradigm(i),:))
end

x = mean(PID(a:z,:)); x = x(:);
y = fA_gain(:);

xx = NaN(length(unique(paradigm)),1);
yy = NaN(length(unique(paradigm)),1);
ww = NaN(length(unique(paradigm)),1);

for i = 1:length(yy)
	xx(i) = mean(x(paradigm==i));
	yy(i) = mean(y(paradigm==i));
	ww(i) = 1./sem(y(paradigm==i));
end


ff = fit(xx(:),yy(:),'power1','Weights',ww);
clear l
l(1) = plot(axes_handles(3),sort(mean(PID(a:z,:))),ff(sort(mean(PID(a:z,:)))),'k--');
L = {};
L{1} = ['y = \alpha x^{', oval(ff.b), '}'];


fo = fitoptions('power1');
fo.Upper = [NaN -1];
fo.Lower = [NaN -1];
fo.Weights = ww;
ff = fit(xx(:),yy(:),'power1',fo);
l(2) = plot(axes_handles(3),sort(mean(PID(a:z,:))),ff(sort(mean(PID(a:z,:)))),'k');
L{2} = ['y = \alpha x^{', oval(ff.b), '}'];

legend(l,L)
legend('boxoff')
set(axes_handles(3),'XScale','log','YScale','log')
xlabel(axes_handles(3),'Mean Stimulus (V)')
ylabel(axes_handles(3),'ORN Gain (Hz/V)')

%    ########    ###     ######  ########     ######      ###    #### ##    ## 
%    ##         ## ##   ##    ##    ##       ##    ##    ## ##    ##  ###   ## 
%    ##        ##   ##  ##          ##       ##         ##   ##   ##  ####  ## 
%    ######   ##     ##  ######     ##       ##   #### ##     ##  ##  ## ## ## 
%    ##       #########       ##    ##       ##    ##  #########  ##  ##  #### 
%    ##       ##     ## ##    ##    ##       ##    ##  ##     ##  ##  ##   ### 
%    ##       ##     ##  ######     ##        ######   ##     ## #### ##    ## 

%     ######   #######  ##    ## ######## ########   #######  ##       
%    ##    ## ##     ## ###   ##    ##    ##     ## ##     ## ##       
%    ##       ##     ## ####  ##    ##    ##     ## ##     ## ##       
%    ##       ##     ## ## ## ##    ##    ########  ##     ## ##       
%    ##       ##     ## ##  ####    ##    ##   ##   ##     ## ##       
%    ##    ## ##     ## ##   ###    ##    ##    ##  ##     ## ##       
%     ######   #######  ##    ##    ##    ##     ##  #######  ######## 

clearvars -except axes_handles being_published

p = '/local-data/DA-paper/large-variance-flicker/LFP/';
[PID, LFP, fA, paradigm, orn, fly, AllControlParadigms, paradigm_hashes]  = consolidateData(p,1);


% set to NaN firing rates that are 0
fA(:,max(fA) == 0) = NaN;

% throw out trials where we didn't record the LFP, for whatever reason
not_LFP = 0*orn;
for i = 1:width(LFP)
	not_LFP(i) = abs(mean2(LFP(:,i)));
end
LFP(:,not_LFP< 0.5) = NaN;

% filter the LFP
filteredLFP = LFP;
for i = 1:width(LFP)
	a = find(~isnan(LFP(:,i)),1,'first');
	z = find(~isnan(LFP(:,i)),1,'last');
	if isempty(a)
		a = 1;
	end
	if isempty(z)
		z = length(LFP);
	end
	try
		filteredLFP(a:z,i) = 10*bandPass(LFP(a:z,i),1000,10);
	catch
	end
end

% extract LFP filters
if ~exist('K','var')
	K = NaN(1e3,length(unique(orn)));
	for i = 1:length(unique(orn))
		resp = mean2(filteredLFP(20e3:55e3,orn==i));
		stim = mean2(PID(20e3:55e3,orn==i));
		temp = fitFilter2Data(stim,resp,'reg',1,'filter_length',1400,'offset',400);
		% throw out 200ms on either end
		temp(1:200) = [];
		temp(end-199:end) = [];
		K(:,i) = temp;
	end
end

% extract fA filters
if ~exist('K3','var')
	K3 = NaN(1e3,length(unique(orn)));
	for i = 1:length(unique(orn))
		stim = mean2(PID(20e3:55e3,orn==i));
		resp = mean2(fA(20e3:55e3,orn==i));
		temp = fitFilter2Data(stim,resp,'reg',1,'filter_length',1400,'offset',400);
		% throw out 200ms on either end
		temp(1:200) = [];
		temp(end-199:end) = [];
		K3(:,i) = temp;
	end
end

% make LN predictions of the LFP and the response
time = 1e-3*(1:length(PID));
fp = NaN(length(fA),length(unique(orn)));
for i = 1:length(unique(orn))
	filtertime = 1e-3*(1:1e3)-.2;
	fp(:,i) = convolve(time,mean2(PID(:,orn==i)),K3(:,i),filtertime);
	% correct for some trivial scaling
	a = fp(20e3:55e3,i);
	b = mean2(fA(20e3:55e3,orn == i));
	temp  = isnan(a) | isnan(b);
	temp = fit(a(~temp),b(~temp),'poly5'); 
	fp(:,i) = temp(fp(:,i));
end

LFP_pred = NaN(length(fA),length(unique(orn)));
for i = 1:length(unique(orn))
	filtertime = 1e-3*(1:1e3)-.2;
	LFP_pred(:,i) = convolve(time,mean2(PID(:,orn==i)),K(:,i),filtertime);
	% correct for some trivial scaling
	a = LFP_pred(:,i);
	b = mean2(filteredLFP(:,orn == i));
	temp  = isnan(a) | isnan(b);
	temp = fit(a(~temp),b(~temp),'poly7'); 
	LFP_pred(:,i) = temp(LFP_pred(:,i));
end

history_lengths = logspace(-1,1,30);

for i = 1:length(unique(orn))
	% first we do the firing rates
	clear ph
	ph(3) = axes_handles(6);
	ph(4) = axes_handles(7);
	
	resp = mean2(fA(:,orn==i));
	stim = mean2(PID(:,orn==i));
	pred = (fp(:,i));

	stim = stim(20e3:55e3);
	pred = pred(20e3:55e3);
	resp = resp(20e3:55e3);

	[p,~,~,~,~,history_lengths,handles]=gainAnalysisWrapper('response',resp,'prediction',pred,'stimulus',stim,'time',1e-3*(1:length(resp)),'ph',ph,'history_lengths',history_lengths,'example_history_length',history_lengths(10),'use_cache',1,'engine',@gainAnalysis);
	

	% cosmetics
	h=get(ph(3),'Children');
	rm_this = [];
	for j = 1:length(h)
		if strcmp(get(h(j),'Marker'),'.')
			rm_this = [rm_this j];
		end
	end
	delete(h(rm_this))
	delete(handles.green_line)
	delete(handles.red_line)
	delete(handles.vert_line)

	set(handles.red_dots,'SizeData',516)
	set(handles.green_dots,'SizeData',516)


	% now do the same for the LFP
	clear ph
	ph(3) = axes_handles(4);
	ph(4) = axes_handles(5);

	resp = mean2(filteredLFP(:,orn==i));
	pred = LFP_pred(:,i);
	pred = pred(20e3:55e3);
	resp = resp(20e3:55e3);

	[p,~,~,~,~,history_lengths,handles] = gainAnalysisWrapper('response',resp,'prediction',pred,'stimulus',stim,'time',1e-3*(1:length(resp)),'ph',ph,'history_lengths',history_lengths,'example_history_length',history_lengths(10),'use_cache',1,'engine',@gainAnalysis);
	
	% cosmetics
	h=get(ph(3),'Children');
	rm_this = [];
	for j = 1:length(h)
		if strcmp(get(h(j),'Marker'),'.')
			rm_this = [rm_this j];
		end
	end
	delete(h(rm_this))
	delete(handles.green_line)
	delete(handles.red_line)
	delete(handles.vert_line)

	set(handles.red_dots,'SizeData',516)
	set(handles.green_dots,'SizeData',516)
end


set(axes_handles(5),'XLim',[.1 10])
set(axes_handles(7),'XLim',[.1 10])

title(axes_handles(4),'')
ylabel(axes_handles(4),'\DeltaLFP (mV)')
xlabel(axes_handles(4),'LN prediction (mV)')

title(axes_handles(6),'')
ylabel(axes_handles(6),'Firing Rate (Hz)')
xlabel(axes_handles(6),'LN prediction (Hz)')

prettyFig('fs=15;','lw=1.5;','plw=2;')

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
