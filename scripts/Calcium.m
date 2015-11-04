% Calcium.m
% 
% created by Srinivas Gorur-Shandilya at 4:00 , 04 November 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

%% Calcium Manipulations of Extracellular Medium
% In this document we manipulate the extra-cellular calcium concentration using EDTA (1mM, to lower it), or increasing the Calcium concentration to 10mM (from a normal conc. of 1mM). 


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

[alldata(1).PID, alldata(1).LFP, alldata(1).fA, alldata(1).paradigm, alldata(1).orn, alldata(1).fly, alldata(1).AllControlParadigms, alldata(1).paradigm_hashes] = consolidateData('/local-data/DA-paper/calcium/low/',true);

[alldata(2).PID, alldata(2).LFP, alldata(2).fA, alldata(2).paradigm, alldata(2).orn, alldata(2).fly, alldata(2).AllControlParadigms, alldata(2).paradigm_hashes] = consolidateData('/local-data/DA-paper/calcium/normal/',true);

[alldata(3).PID, alldata(3).LFP, alldata(3).fA, alldata(3).paradigm, alldata(3).orn, alldata(3).fly, alldata(3).AllControlParadigms, alldata(3).paradigm_hashes] = consolidateData('/local-data/DA-paper/calcium/high/',true);

for ai = 1:length(alldata)
	% remove baseline from all PIDs
	for i = 1:width(alldata(ai).PID)
		alldata(ai).PID(:,i) = alldata(ai).PID(:,i) - mean(alldata(ai).PID(1:5e3,i));
	end

	% remove baseline from all LFPs
	for i = 1:width(alldata(ai).LFP)
		alldata(ai).LFP(:,i) = alldata(ai).LFP(:,i) - mean(alldata(ai).LFP(1:5e3,i));
	end

	% band pass all the LFP
	alldata(ai).filtered_LFP = alldata(ai).LFP;
	for i = 1:width(alldata(ai).LFP)
		alldata(ai).filtered_LFP(:,i) = bandPass(alldata(ai).LFP(:,i),1000,10);
	end

	% remove "Flicker" from paradigm names
	for i = 1:length(alldata(ai).AllControlParadigms)
		alldata(ai).AllControlParadigms(i).Name = strrep(alldata(ai).AllControlParadigms(i).Name,'Flicker-','');
	end

end

%% Baseline Firing
% In this section we look at the baseline firing of the neuron. 

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
for ai = 1:length(alldata)
	temp = nonnans(nonzeros(mean(alldata(ai).fA(1:5e3,:))));
	errorbar(ai,mean(temp),sem(temp),'k')
end
set(gca,'XTick',1:length(alldata),'XTickLabel',{'Low Ca^{2+}','Normal Ca^{2+}','High Ca^{2+}'},'XTickLabelRotation',45)
ylabel('Baseline Firing Rate (Hz)')
prettyFig()

if being_published
	snapnow
	delete(gcf)
end


%% LFP Gain
% In this section we compare the LFP response to a fluctuating stimulus between normal calcium and low and high manipulations. 

a = 10e3;
z = 45e3;
ss = 20;
figure('outerposition',[0 0 1200 700],'PaperUnits','points','PaperSize',[1200 700]); hold on
for i = 1:3
	subplot(2,3,i), hold on
	normal_LFP =  alldata(2).filtered_LFP(a:z,alldata(2).paradigm==i);
	test_LFP = alldata(1).filtered_LFP(a:z,alldata(1).paradigm==i);
	m = 1.2*min([test_LFP(:); normal_LFP(:)]);
	M = 1.2*max([test_LFP(:); normal_LFP(:)]);
	plot([m M],[m M],'k--')
	plot(normal_LFP(1:ss:end),test_LFP(1:ss:end),'b.')
	axis square
	xlabel('\Delta LFP at normal Ca^{2+}')
	ylabel('\Delta LFP at low Ca^{2+}')
	title(alldata(1).AllControlParadigms(find(alldata(2).paradigm==i)).Name)

	subplot(2,3,3+i), hold on
	test_LFP = alldata(3).filtered_LFP(a:z,alldata(3).paradigm==i);
	m = 1.2*min([test_LFP(:); normal_LFP(:)]);
	M = 1.2*max([test_LFP(:); normal_LFP(:)]);
	plot([m M],[m M],'k--')
	plot(normal_LFP(1:ss:end),test_LFP(1:ss:end),'r.')
	axis square
	xlabel('\Delta LFP at normal Ca^{2+}')
	ylabel('\Delta LFP at high Ca^{2+}')
end

prettyFig()

if being_published
	snapnow
	delete(gcf)
end

%% Firing Gain 
% In this section we compare the Firing response to a fluctuating stimulus between normal calcium and low and high manipulations. 

a = 10e3;
z = 45e3;
ss = 20;
figure('outerposition',[0 0 1200 700],'PaperUnits','points','PaperSize',[1200 700]); hold on
for i = 1:3
	subplot(2,3,i), hold on
	normal_LFP =  alldata(2).fA(a:z,alldata(2).paradigm==i);
	test_LFP = alldata(1).fA(a:z,alldata(1).paradigm==i);
	m = 0;
	M = 100;
	plot([m M],[m M],'k--')
	plot(normal_LFP(1:ss:end),test_LFP(1:ss:end),'b.')
	axis square
	xlabel('Firing Rate at normal Ca^{2+}')
	ylabel('Firing Rate at low Ca^{2+}')
	title(alldata(1).AllControlParadigms(find(alldata(2).paradigm==i)).Name)

	subplot(2,3,3+i), hold on
	test_LFP = alldata(3).fA(a:z,alldata(3).paradigm==i);
	m = 0;
	M = 100;
	plot([m M],[m M],'k--')
	plot(normal_LFP(1:ss:end),test_LFP(1:ss:end),'r.')
	axis square
	xlabel('Firing Rate at normal Ca^{2+}')
	ylabel('Firing Rate at high Ca^{2+}')
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

t = toc;

%% 
% This document was built in: 
disp(strcat(oval(t,3),' seconds.'))

% tag the file as being published 

if being_published
	unix(['tag -a published ',which(mfilename)]);
	unix(['tag -r publish-failed ',which(mfilename)]);
end
