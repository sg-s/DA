% BaselineFiring.m
% Analysis of firing properties of ORNs at rest
% 
% created by Srinivas Gorur-Shandilya at 3:51 , 30 July 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


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

%% Baseline Firing 
% In this document, we carefully investigate the firing statistics of the A and B neurons at rest, i.e., with no stimulation whatsoever. We recorded from ab2 and ab3 sensilla. In this entire document, A is plotted in red, and B is plotted in blue. 

A = [];
B = [];
orn = [];
pathname = '/local-data/orn/baseline-firing/ab3/';
allfiles = dir([pathname '*.mat']);
for i = 1:length(allfiles)
	load(strcat(pathname,allfiles(i).name));
	orn = [orn; i*ones(width(spikes.A),1)];
	A = [A  spikes.A'];
	B = [B  spikes.B'];
end
ab3.A = A;
ab3.B = B;
ab3.orn = orn;


A = [];
B = [];
orn = [];
pathname = '/local-data/orn/baseline-firing/ab2/';
allfiles = dir([pathname '*.mat']);
for i = 1:length(allfiles)
	load(strcat(pathname,allfiles(i).name));
	A = [A  spikes.A'];
	orn = [orn; i*ones(width(spikes.A),1)];
	B = [B  spikes.B'];
end
ab2.A = A;
ab2.B = B;
ab2.orn = orn;

%    ###    ########   #######  
%   ## ##   ##     ## ##     ## 
%  ##   ##  ##     ##        ## 
% ##     ## ########   #######  
% ######### ##     ## ##        
% ##     ## ##     ## ##        
% ##     ## ########  ######### 


%% ab2 
% First we analyse the firing properties of the ab2 sensilla at rest. The following figure shows the raster of all the spikes we are going to analyse. (A is in red, B is in blue.)

figure('outerposition',[0 0 1500 700],'PaperUnits','points','PaperSize',[1500 700]); hold on
raster2(A,B)
xlabel('Time (s)')
title('ab2')
prettyFig('plw=.5;')

if being_published
	snapnow
	delete(gcf)
end

% In the following figure, we plot the inter-spike interval (ISI) histogram for the A and the B neurons. The different lines indicate different neurons, and the shading is the s.e.m. (A is in red, B is in blue.)


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
for i = 1:max(orn)
	pt = orn == i;
	errorShade(hx,mean2(isiB(:,pt)),sem(isiB(:,pt)),'Color',[0 0 1]);
	errorShade(hx,mean2(isiA(:,pt)),sem(isiA(:,pt)),'Color',[1 0 0]);
end
xlabel('ISI (ms)')
ylabel('Probability')
title('ab2')
prettyFig()

if being_published
	snapnow
	delete(gcf)
end


%%
% Now we plot the probability of B spiking conditional on a A spike, and vice versa. This is a measure of how much the two neurons affect each other. Note that this is not the same as the ISI distribution, because this depends on all spikes, not just adjacent pairs of spikes. In colours are the self-probability curves (e.g., prob. of A spike given A spike), and in black are the cross-probability curves (e.g., prob of A spike given B spike). 

pAA = zeros(100,width(A));
pBB = zeros(100,width(A));
pAB = zeros(100,width(A));
pBA = zeros(100,width(A));
for i = 1:width(A)
	[pAA(:,i),x] = condSpikeProb(A(:,i),A(:,i),1,1e-2);
	[pBB(:,i),x] = condSpikeProb(B(:,i),B(:,i),1,1e-2);
	[pAB(:,i),x] = condSpikeProb(A(:,i),B(:,i),1,1e-2);
	[pBA(:,i),x] = condSpikeProb(B(:,i),A(:,i),1,1e-2);
end

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
clear l

temp = shadedErrorBar(x,mean2(pBB),sem(pBB),{'Color',[0 0 1]});
l(1) = temp.mainLine;
temp = shadedErrorBar(x,mean2(pAB),sem(pAB),{'Color',[0 0 0]});
l(2) = temp.mainLine;

legend(l,{'conditional on B spike','conditional on A spike'})
xlabel('Time since spike (s)')
ylabel('p(B spike observed)')
clear l
subplot(1,2,2), hold on

temp = shadedErrorBar(x,mean2(pAA),sem(pAA),{'Color',[1 0 0]});
l(1) = temp.mainLine;
temp = shadedErrorBar(x,mean2(pBA),sem(pBA),{'Color',[0 0 0]});
l(2) = temp.mainLine;

legend(l,{'conditional on A spike','conditional on B spike'})
xlabel('Time since spike (s)')
ylabel('p(A spike observed)')
suptitle('ab2 Sensilla')
prettyFig()

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
% Now we look at ab3. The following figure shows the raster of all the spikes we are going to analyse. (A is in red, B is in blue.)

A = ab3.A; B = ab3.B; orn = ab3.orn;

figure('outerposition',[0 0 1500 700],'PaperUnits','points','PaperSize',[1500 700]); hold on
raster2(A,B)
xlabel('Time (s)')
title('ab3')
prettyFig('plw=.5;')

if being_published
	snapnow
	delete(gcf)
end

% In the following figure, we plot the inter-spike interval (ISI) histogram for the A and the B neurons. The different lines indicate different neurons, and the shading is the s.e.m. (A is in red, B is in blue.)


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
for i = 1:max(orn)
	pt = orn == i;
	errorShade(hx,mean2(isiB(:,pt)),sem(isiB(:,pt)),'Color',[0 0 1]);
	errorShade(hx,mean2(isiA(:,pt)),sem(isiA(:,pt)),'Color',[1 0 0]);
end
xlabel('ISI (ms)')
ylabel('Probability')
title('ab3')
prettyFig()

if being_published
	snapnow
	delete(gcf)
end
%%
% Now we plot the probability of B spiking conditional on a A spike, and vice versa. This is a measure of how much the two neurons affect each other. Note that this is not the same as the ISI distribution, because this depends on all spikes, not just adjacent pairs of spikes. In colours are the self-probability curves (e.g., prob. of A spike given A spike), and in black are the cross-probability curves (e.g., prob of A spike given B spike). 

pAA = zeros(100,width(A));
pBB = zeros(100,width(A));
pAB = zeros(100,width(A));
pBA = zeros(100,width(A));
for i = 1:width(A)
	[pAA(:,i),x] = condSpikeProb(A(:,i),A(:,i),1,1e-2);
	[pBB(:,i),x] = condSpikeProb(B(:,i),B(:,i),1,1e-2);
	[pAB(:,i),x] = condSpikeProb(A(:,i),B(:,i),1,1e-2);
	[pBA(:,i),x] = condSpikeProb(B(:,i),A(:,i),1,1e-2);
end

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
clear l

temp = shadedErrorBar(x,mean2(pBB),sem(pBB),{'Color',[0 0 1]});
l(1) = temp.mainLine;
temp = shadedErrorBar(x,mean2(pAB),sem(pAB),{'Color',[0 0 0]});
l(2) = temp.mainLine;

legend(l,{'conditional on B spike','conditional on A spike'},'location','southeast')
xlabel('Time since spike (s)')
ylabel('p(B spike observed)')
clear l
subplot(1,2,2), hold on

temp = shadedErrorBar(x,mean2(pAA),sem(pAA),{'Color',[1 0 0]});
l(1) = temp.mainLine;
temp = shadedErrorBar(x,mean2(pBA),sem(pBA),{'Color',[0 0 0]});
l(2) = temp.mainLine;

legend(l,{'conditional on A spike','conditional on B spike'},'location','southeast')
xlabel('Time since spike (s)')
ylabel('p(A spike observed)')
suptitle('ab3 Sensilla')
prettyFig()

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
% add homebrew path
path1 = getenv('PATH');
path1 = [path1 ':/usr/local/bin'];
setenv('PATH', path1);

if being_published
	unix(['tag -a published ',which(mfilename)]);
	unix(['tag -r publish-failed ',which(mfilename)]);
end

