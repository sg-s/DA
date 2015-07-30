% BaselineFiring.m
% 
% created by Srinivas Gorur-Shandilya at 3:51 , 30 July 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


% this code determines if this function is being called by publish() or not
calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end
tic

%% Baseline Firing 
% In this document, we carefully investigate the firing statistics of the A and B neurons at rest, i.e., with no stimulation whatsoever. We recorded from ab2 and ab3 sensilla. 

A = [];
B = [];
pathname = '/local-data/orn/baseline-firing/ab3/';
allfiles = dir([pathname '*.mat']);
for i = 1:length(allfiles)
	load(strcat(pathname,allfiles(i).name));
	A = [A  spikes.A'];
	B = [B  spikes.B'];
end
ab3.A = A;
ab3.B = B;


A = [];
B = [];
pathname = '/local-data/orn/baseline-firing/ab2/';
allfiles = dir([pathname '*.mat']);
for i = 1:length(allfiles)
	load(strcat(pathname,allfiles(i).name));
	A = [A  spikes.A'];
	B = [B  spikes.B'];
end
ab2.A = A;
ab2.B = B;

%    ###    ########   #######  
%   ## ##   ##     ## ##     ## 
%  ##   ##  ##     ##        ## 
% ##     ## ########   #######  
% ######### ##     ## ##        
% ##     ## ##     ## ##        
% ##     ## ########  ######### 


%% ab2 
% First we analyse the firing properties of the ab2 sensilla at rest. The following figure shows the raster of all the spikes we are going to analyse:

figure('outerposition',[0 0 1500 700],'PaperUnits','points','PaperSize',[1500 700]); hold on
raster2(A,B)
xlabel('Time (s)')
title('ab2')
PrettyFig('plw=.5;')

if being_published
	snapnow
	delete(gcf)
end

% In the following figure, we plot the inter-spike interval (ISI) histogram for the A and the B neurons: 



figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
hx = 0:10:3e3;
isiA = zeros(length(hx),width(A));
for i = 1:width(A)
	y = diff(find(A(:,i)))/10; % in ms
	hy = hist(y,hx); % up to 1 s
	hy = hy/sum(hy);
	isiA(:,i) = hy;
end
isiB = zeros(length(hx),width(B));
for i = 1:width(B)
	y = diff(find(B(:,i)))/10; % in ms
	hy = hist(y,hx); % up to 1 s
	hy = hy/sum(hy);
	isiB(:,i) = hy;
end
errorShade(hx,mean2(isiB),sem(isiB),'Color',[0 0 1]);
errorShade(hx,mean2(isiA),sem(isiA),'Color',[1 0 0]);
xlabel('ISI (ms)')
ylabel('Probability')
title('ab2')
PrettyFig()

if being_published
	snapnow
	delete(gcf)
end

%                             ###    ########   #######  
%                            ## ##   ##     ## ##     ## 
%                           ##   ##  ##     ##        ## 
%                          ##     ## ########   #######  
%                          ######### ##     ##        ## 
%                          ##     ## ##     ## ##     ## 
%                          ##     ## ########   #######  



%% ab3
% Now we look at ab3. The following figure shows the raster of all the spikes we are going to analyse:

A = ab3.A; B = ab3.B;

figure('outerposition',[0 0 1500 700],'PaperUnits','points','PaperSize',[1500 700]); hold on
raster2(A,B)
xlabel('Time (s)')
title('ab3')
PrettyFig('plw=.5;')

if being_published
	snapnow
	delete(gcf)
end

% In the following figure, we plot the inter-spike interval (ISI) histogram for the A and the B neurons: 


figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
hx = 0:10:3e3;
isiA = zeros(length(hx),width(A));
for i = 1:width(A)
	y = diff(find(A(:,i)))/10; % in ms
	hy = hist(y,hx); % up to 1 s
	hy = hy/sum(hy);
	isiA(:,i) = hy;
end
isiB = zeros(length(hx),width(B));
for i = 1:width(B)
	y = diff(find(B(:,i)))/10; % in ms
	hy = hist(y,hx); % up to 1 s
	hy = hy/sum(hy);
	isiB(:,i) = hy;
end
errorShade(hx,mean2(isiB),sem(isiB),'Color',[0 0 1]);
errorShade(hx,mean2(isiA),sem(isiA),'Color',[1 0 0]);
xlabel('ISI (ms)')
ylabel('Probability')
title('ab3')
PrettyFig()

if being_published
	snapnow
	delete(gcf)
end

%%
% So there is something strikingly different about the B neuron in ab2 and ab3. In ab3, it looks normal (Poisson-distributed), just like the A neuron, but in the ab3 sensilla, it has a weird hump at 700ms+. 


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

% tag the file as being published 
% add homebrew path
path1 = getenv('PATH');
path1 = [path1 ':/usr/local/bin'];
setenv('PATH', path1);

if being_published
	unix(strjoin({'tag -a published',which(mfilename)}));
end
