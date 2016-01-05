% obp28a_analysis_2.m
% 
% created by Srinivas Gorur-Shandilya at 3:42 , 01 January 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.ts


%% OBP28a analysis of Gain
% In this document, we see if the gain properties of the ab8A ORN are changed when we delete the OBP 28a. 



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

% prep and merge data
use_cache = 1;
[PID, LFP, fA, paradigm, orn,fly] = consolidateData('/local-data/obp/ab8/ko',use_cache);
geno = ones(length(orn),1);

[PID2, LFP2, fA2, paradigm2, orn2,fly2] = consolidateData('/local-data/obp/ab8/wcs',use_cache);
geno2 = 2*ones(length(orn2),1);

PID = [PID PID2]; clear PID2
fA = [fA fA2]; clear fA2
LFP = [LFP LFP2]; clear LFP2
paradigm = [paradigm paradigm2]; clear paradigm2
orn = [orn orn2]; clear orn2
fly = [fly; fly2]; clear fly2
geno = [geno; geno2]; clear geno2
paradigm = paradigm(:);
orn = orn(:);


% clean up data
rm_this = isnan(sum(LFP));
fA(:,rm_this) = [];
PID(:,rm_this) = [];
LFP(:,rm_this) = [];
paradigm(rm_this) = [];
orn(rm_this) = [];
fly(rm_this) = [];
geno(rm_this) = [];
fA(:,sum(fA) == 0) = NaN;


% bandPass LFP
filtered_LFP = LFP;
for i = 1:width(LFP)
	filtered_LFP(:,i) = bandPass(LFP(:,i),1000,10);
end

a = 10e3; z = 50e3;
[~,~,LFP_gain] = extractFilters(PID,filtered_LFP,'use_cache',true,'a',a,'z',z);

[~,~,fA_gain] = extractFilters(PID,fA,'use_cache',true,'a',a,'z',z);

%% Comparison of Gain 
% In the following figure we compare the gain in the LFP and in the firing rate between the wCS control (black) and the CRISPR OBP deletion (red). 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
x = mean(PID(a:z,geno==1)); y = LFP_gain(geno==1);
ff = fit(x(:),y(:),'power1');
plot(.1:0.1:2,ff(.1:0.1:2),'r')
plot(x,y,'r+');

x = mean(PID(a:z,geno==2)); y = LFP_gain(geno==2);
ff = fit(x(:),y(:),'power1');
plot(.1:0.1:2,ff(.1:0.1:2),'k')
plot(x,y,'k+');

set(gca,'XScale','log','YScale','log')
xlabel('Mean Stimulus (V)')
ylabel('LFP Gain (mV/V)')

subplot(1,2,2), hold on
x = mean(PID(a:z,geno==1)); y = fA_gain(geno==1);
x(isnan(y)) = []; y(isnan(y)) = [];
ff = fit(x(:),y(:),'power1');
plot(.1:0.1:2,ff(.1:0.1:2),'r')
plot(x,y,'r+');

x = mean(PID(a:z,geno==2)); y = fA_gain(geno==2);
x(isnan(y)) = []; y(isnan(y)) = [];
ff = fit(x(:),y(:),'power1');
plot(.1:0.1:2,ff(.1:0.1:2),'k')
plot(x,y,'k+');

set(gca,'XScale','log','YScale','log')
xlabel('Mean Stimulus (V)')
ylabel('Firing Gain (Hz/V)')


prettyFig('fs=20;')

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

