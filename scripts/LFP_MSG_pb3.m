% Mean Shifted Gaussians + Local Field Potential Measurements in pb3
% 
% created by Srinivas Gorur-Shandilya at 10:40 , 28 October 2015. Contact me at http://srinivas.gs/contact/
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


%           ########  ########       #######  
%           ##     ## ##     ##     ##     ## 
%           ##     ## ##     ##            ## 
%           ########  ########       #######  
%           ##        ##     ##            ## 
%           ##        ##     ##     ##     ## 
%           ##        ########       #######  


%% Gain changes in LFP in pb3
% In this document we measure LFPs from the pb3 sensillum as we excite it with isoamyl acetate. We want to know if the gain scaling exponent in the LFP is still much steeper than -1 (as we saw in ab3). If this is the case, we are more sure that this is real, and not an artifact of dense sensillar interference in the antenna on the LFP. 

[pb3_PID, pb3_LFP, ~, pb3_paradigm, pb3_orn, pb3_fly, pb3_AllControlParadigms] = consolidateData('/local-data/DA-paper/palp/pb3/',1);


% remove baseline from all PIDs
for i = 1:width(pb3_PID)
	pb3_PID(:,i) = pb3_PID(:,i) - mean(pb3_PID(1:5e3,i));
end

% sort the paradigms sensibly
sort_value = [];
for i = 1:length(pb3_AllControlParadigms)
	sort_value(i) = (mean(pb3_AllControlParadigms(i).Outputs(1,:)));
end

[~,idx] = sort(sort_value);
pb3_AllControlParadigms = pb3_AllControlParadigms(idx);
paradigm_new = pb3_paradigm*NaN;
for i = 1:length(idx)
	paradigm_new(pb3_paradigm == idx(i)) = i;
end
pb3_paradigm = paradigm_new;

% remove "Flicker" from paradigm names
for i = 1:length(pb3_AllControlParadigms)
	pb3_AllControlParadigms(i).Name = strrep(pb3_AllControlParadigms(i).Name,'Flicker-','');
end

% throw out trials where we didn't record the LFP, for whatever reason
not_LFP = find((max(abs(pb3_LFP))) < 0.1);
pb3_LFP(:,not_LFP) = NaN;

% throw our bad traces
bad_trials =  (isnan(sum(pb3_LFP)));
pb3_LFP(:,bad_trials) = [];
pb3_PID(:,bad_trials) = [];
pb3_paradigm(bad_trials) = [];
pb3_orn(bad_trials) = [];

% band pass all the LFP
for i = 1:width(pb3_LFP)
	pb3_LFP(:,i) = bandPass(pb3_LFP(:,i),1000,10);
end

a = 10e3; z = 50e3;
[pb3_K,pb3_LFP_pred,pb3_LFP_gain,pb3_LFP_gain_err] = extractFilters(pb3_PID,pb3_LFP,'use_cache',true,'a',a,'z',z);

figure('outerposition',[0 0 800 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(nanmean(pb3_PID(a:z,:)),pb3_LFP_gain,'k+')

x = mean(pb3_PID(a:z,:)); x = x(:);
y = pb3_LFP_gain(:);

xx = NaN(length(unique(pb3_paradigm)),1);
yy = NaN(length(unique(pb3_paradigm)),1);
ww = NaN(length(unique(pb3_paradigm)),1);

for i = 1:length(yy)
	xx(i) = mean(x(pb3_paradigm==i));
	yy(i) = mean(y(pb3_paradigm==i));
	ww(i) = 1./sem(y(pb3_paradigm==i));
end
ww(isinf(ww)) = max(ww(~isinf(ww)));
ff = fit(xx(:),yy(:),'power1','Weights',ww);
clear l
l(1) = plot(sort(xx),ff(sort(xx)),'k--');
L = {};
L{1} = ['y = \alpha (x^\beta) ,\beta = ', oval(ff.b) ,' , r^2=' oval(rsquare(ff(x),y))];

fo = fitoptions('power1');
fo.Upper = [NaN -1];
fo.Lower = [NaN -1];
fo.Weights = ww;
ff = fit(xx(:),yy(:),'power1',fo);
l(2) = plot(sort(xx),ff(sort(xx)),'k');
L{2} = ['y = \alpha (x^\beta) ,\beta := -1, r^2=' oval(rsquare(ff(x),y))];

legend(l,L,'Location','eastoutside')
set(gca,'XScale','log','YScale','log','XLim',[.05 .3],'YLim',[min(pb3_LFP_gain)/2 max(pb3_LFP_gain)*2])
xlabel('Mean Stimulus (V)')
ylabel('LFP Gain (mV/V)')

prettyFig;

if being_published
	snapnow
	delete(gcf)
end

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

%%
% This file has the following external dependencies:
showDependencyHash(mfilename);

t = toc;

%% 
% This document was built in: 
disp(strcat(oval(t,3),' seconds.'))

% tag the file as being published 

if being_published
	unix(['tag -a published ',which(mfilename)]);
	unix(['tag -r publish-failed ',which(mfilename)]);
end


