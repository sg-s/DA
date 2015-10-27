% GAlpha.m
% 
% created by Srinivas Gorur-Shandilya at 8:04 , 26 October 2015. Contact me at http://srinivas.gs/contact/
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


%% Stimulus
% In this section we show the stimulus. 

[PID, LFP, fA, paradigm, orn, fly, AllControlParadigms, paradigm_hashes] = consolidateData('/local-data/DA-paper/g-alpha/rnai',1);



% remove baseline from all PIDs
for i = 1:width(PID)
	PID(:,i) = PID(:,i) - mean(PID(1:5e3,i));
end

% remove baseline from all LFPs
for i = 1:width(LFP)
	LFP(:,i) = LFP(:,i) - mean(LFP(1:5e3,i));
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

% throw out trials where i think we lost the neuron
lost_neuron = ones(width(LFP),1);
for i = 1:width(LFP)
	temp = LFP(50e3:60e3,i);
	lost_neuron(i) = nanstd(temp);
end
LFP(:,lost_neuron<.1) = NaN;


% throw our bad traces
bad_trials = isnan(sum(LFP));
LFP(:,bad_trials) = [];
PID(:,bad_trials) = [];
fA(:,bad_trials) = [];
paradigm(bad_trials) = [];
orn(bad_trials) = [];

% band pass all the LFP
filtered_LFP = LFP;
for i = 1:width(LFP)
	filtered_LFP(:,i) = bandPass(LFP(:,i),1000,10);
end


figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,4,1:3), hold on
c = parula(1+length(unique(paradigm)));
for i = 1:length(unique(paradigm))
	plot_this = PID(40e3:45e3,paradigm==i);
	time = 40+1e-3*(1:length(plot_this));
	plot(time,mean(plot_this,2),'Color',c(i,:),'LineWidth',2);
end
ylabel('Stimulus (V)')
xlabel('Time (s)')
set(gca,'YLim',[0 2.5])

subplot(1,4,4), hold on
c = parula(1+length(unique(paradigm)));
for i = 1:length(unique(paradigm))
	hist_this = PID(20e3:55e3,paradigm==i);
	xx =  linspace(min(min(hist_this)),max(max(hist_this)),50);
	y = NaN(sum(paradigm==i),50);
	for j = 1:sum(paradigm==i)
		y(j,:) = hist(hist_this(:,j),xx);
		y(j,:) = y(j,:)/sum(y(j,:));
	end
	if width(y)> 1
		plot(mean2(y),(xx),'Color',c(i,:));
	else
		plot(y,xx,'Color',c(i,:))
	end
end

xlabel('p(stimulus)')
set(gca,'YLim',[0 2.5])
prettyFig('fs=14;');

if being_published
	snapnow
	delete(gcf)
end

%% LFP Gain
% In this section we analyse the gain in the LFP


a = 10e3; z = 50e3;
[K1,LFP_pred,LFP_gain,LFP_gain_err] = extractFilters(PID,filtered_LFP,'use_cache',true,'a',a,'z',z);


ss  =50;
c= parula(max(paradigm)+1);
l = [];
filtertime = 1e-3*(1:length(K1))-.1;
figure('outerposition',[0 0 1400 500],'PaperUnits','points','PaperSize',[1400 500]); hold on
subplot(1,10,1:5), hold on
for i = 1:max(paradigm)
	plot_this = find(paradigm == i);
	plot_this = setdiff(plot_this,find(isnan(sum(K1))));
	y = filtered_LFP(a:z,plot_this);
	if width(y) > 1
		y = mean2(y);
		x = mean2(LFP_pred(a:z,plot_this));
	end
	
	l(i) = plot(x(1:ss:end),y(1:ss:end),'.','Color',c(i,:));
end
legend(l,{AllControlParadigms.Name},'Location','southeastoutside')
xlabel('Linear Prediction')
ylabel('\DeltaLFP (mV)')

subplot(1,10,7:10), hold on
for i = 1:width(LFP)
	plot(mean(PID(a:z,i)),LFP_gain(i),'+','Color',c(paradigm(i),:))
end

% x = mean(PID(a:z,:)); x = x(:);
% y = LFP_gain(:);

% xx = NaN(length(unique(paradigm)),1);
% yy = NaN(length(unique(paradigm)),1);
% ww = NaN(length(unique(paradigm)),1);

% for i = 1:length(yy)
% 	xx(i) = mean(x(paradigm==i));
% 	yy(i) = mean(y(paradigm==i));
% 	ww(i) = 1./sem(y(paradigm==i));
% end
% ww(isinf(ww)) = max(ww(~isinf(ww)));

% ff = fit(xx(:),yy(:),'power1','Weights',ww);
% clear l
% l(1) = plot(sort(xx),ff(sort(xx)),'k--');
% L = {};
% L{1} = ['y = \alpha (x^\beta) ,\beta = ', oval(ff.b),', r^2=' oval(rsquare(ff(x),y))];


% fo = fitoptions('power1');
% fo.Upper = [NaN -1];
% fo.Lower = [NaN -1];
% fo.Weights = ww;
% ff = fit(xx(:),yy(:),'power1',fo);
% l(2) = plot(sort(xx),ff(sort(xx)),'k');
% L{2} = ['y = \alpha (x^\beta) ,\beta := -1, r^2=' oval(rsquare(ff(x),y))];

% legend(l,L)
set(gca,'XScale','log','YScale','log')
xlabel('Mean Stimulus (V)')
ylabel('LFP Gain (mV/V)')
set(gca,'YLim',[.1 5])

prettyFig;

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

