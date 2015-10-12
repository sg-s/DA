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
	filtered_LFP(:,i) = 100*bandPass(LFP(:,i),1000,10); % now in mV
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
l(1) = plot(axes_handles(2),sort(xx),ff(sort(xx)),'k--');
L = {};
L{1} = ['y = \alpha x^{', oval(ff.b), '}'];

fo = fitoptions('power1');
fo.Upper = [NaN -1];
fo.Lower = [NaN -1];
fo.Weights = ww;
ff = fit(xx(:),yy(:),'power1',fo);
l(2) = plot(axes_handles(2),sort(xx),ff(sort(xx)),'k');
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
l(1) = plot(axes_handles(3),sort(xx),ff(sort(xx)),'k--');
L = {};
L{1} = ['y = \alpha x^{', oval(ff.b), '}'];


fo = fitoptions('power1');
fo.Upper = [NaN -1];
fo.Lower = [NaN -1];
fo.Weights = ww;
ff = fit(xx(:),yy(:),'power1',fo);
l(2) = plot(axes_handles(3),sort(xx),ff(sort(xx)),'k');
L{2} = ['y = \alpha x^{', oval(ff.b), '}'];

legend(l,L)
legend('boxoff')
set(axes_handles(3),'XScale','log','YScale','log')
xlabel(axes_handles(3),'Mean Stimulus (V)')
ylabel(axes_handles(3),'ORN Gain (Hz/V)')


prettyFig('fs=15;')

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
