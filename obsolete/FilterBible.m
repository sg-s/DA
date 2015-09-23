% FilterBible.m
% 
% created by Srinivas Gorur-Shandilya at 10:54 , 24 June 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end
tic


%% Filter Bible.m
% This document rigorously tests different methods of filter extraction with different kinds of synthetic data. 


%% Performance
% In this section, we see how our different methods perform with problem size in a no-noise, perfect case. 



ss = logspace(3,5.2,20);
nrep  = 3;
time1 = NaN(length(ss),nrep);
K = (filter_alpha2(30,70,1,.2,1:300));

for i = 1:length(ss)
	for j = 1:nrep
		x = randn(floor(ss(i)),1);
		y = filter(K,1,x);
		tic
		FitFilter2Data(x,y,[],'reg=0;');
		time1(i,j) = toc;
		tic
		revCorrFilter(x,y,'reg',0);
		time2(i,j) = toc;
	end
end

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
errorbar(ss,mean(time1'),std(time1'),'Color',[1 0 0]);
errorbar(ss,mean(time2'),std(time2'),'Color',[0 0 1]);
xlabel('System Size')
ylabel('Solution Time (s)')
set(gca,'XScale','log','YScale','log')
legend({'Least Squares','reverse Correlation'},'Location','southeast')
PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%% Regularisation
% Performance of filter extraction techniques is interesting only in cases with non-white inputs. How we regularise is key. Here, we employ ridge regression (adding a constant times an identity matrix to the matrix that is being inverted). For least squares, it is the covariance matrix. For reverse correlation, it is the Toeplitz. 

ss = 1e4;
input_corr = linspace(1,100,4);
nrep  = 3;
q1 = NaN(length(input_corr),nrep);
q2 = NaN(length(input_corr),nrep);
K = (filter_alpha2(30,70,1,.2,1:300));

for i = 1:length(input_corr)
	for j = 1:nrep
		x = randn(ss,1);
		x = filter(ones(input_corr(i),1)/input_corr(i),1,x);
		y = filter(K,1,x);
		tic
		FitFilter2Data(x,y,[],'reg=0;');
		time1(i,j) = toc;
		tic
		revCorrFilter(x,y,'reg',0);
		time2(i,j) = toc;
	end
end



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



