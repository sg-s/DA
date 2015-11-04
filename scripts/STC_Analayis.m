% STC_Analayis.m
%
% created by Srinivas Gorur-Shandilya at 11:08 , 02 November 2015. Contact me at http://srinivas.gs/contact/
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


%% STC Analysis
% In this document we analyse the spike-triggered covariance using some synthetic data. Here, we generate stimuli using some random Gaussian noise, use an exponential filter and a threshold to generate spikes, and then look at the STC of that. 

% generate correlated gaussian white noise
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',1984)); 
S = randn(1e4,1);
S = filtfilt(ones(10,1),10,S);

% Compute filter responses -------------
K = filter_exp(20,1,1:500);
K = K/max(K);
R = filter(K,1,S);
R = R/std(R);
R(R<0) = 0; 
R = (R.*rand(length(R),1));
R(R<std(R)) = 0;
R(R>0) = 1;
R = sparse(R);

analyseSTC(S,R,K)

prettyFig()

if being_published
	snapnow
	delete(gcf)
end

%% A more complicated filter
% Now we repeat the analysis with a more realistic filter:

p.tau1 = 20;
p.tau2 = 40;
p.A = 1;
p.n = 2;
K = filter_gamma2(1:500,p);
K = K/max(K);
R = filter(K,1,S);
R = R/std(R);
R(R<0) = 0; 
R = (R.*rand(length(R),1));
R(R<std(R)) = 0;
R(R>0) = 1;
R = sparse(R);
analyseSTC(S,R,K)

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
