% DivisiveOrSubtractive.m
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


% internal housekeeping: determine if being called by publish or not
calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end
% check if there is a data cache and load it
if exist(strcat(mfilename,'.mat'))
	load(strcat(mfilename,'.mat'))
end


%% Divisive of Subtractive Modulation of Gain? 
% Consider a sensor with some input-output curve shown in the following figure:

figure('outerposition',[0 0 400 400],'PaperUnits','points','PaperSize',[1000 400]); hold on
x = logspace(-4,1);
y = hill([100 1e-1 2],x);
plot(x,y,'k')
set(gca,'XScale','log')
xlabel('Stimulus')
ylabel('Response')

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%% 
% when this sensors response properties are modified in some way (for example, when it is "adapted" to something), its input-output curve can change in the following ways:

figure('outerposition',[0 0 1400 500],'PaperUnits','points','PaperSize',[1400 500]); hold on
subplot(1,3,1), hold on
title('Purely subtractive ')
c=  jet(5);
Kd = logspace(-2,-1,5);
for i = 1:length(c)
	y = hill([100 Kd(i) 2],x);
	plot(x,y,'Color',c(i,:))
end
set(gca,'XScale','log','YLim',[0 100])
xlabel('Stimulus')
ylabel('Response')


subplot(1,3,2), hold on
title('Purely divisive ')
c=  jet(5);
n = linspace(2,10,5);
for i = 1:length(c)
	y = hill([100 1e-1 n(i)],x);
	plot(x,y,'Color',c(i,:))
end
set(gca,'XScale','log','YLim',[0 100])
xlabel('Stimulus')
ylabel('Response')


subplot(1,3,3), hold on
title('Mixed')
c=  jet(5);
n = linspace(2,10,5);
Kd = logspace(-2,-1,5);
for i = 1:length(c)
	y = hill([100 Kd(i) n(i)],x);
	plot(x,y,'Color',c(i,:))
end
set(gca,'XScale','log','YLim',[0 100])
xlabel('Stimulus')
ylabel('Response')


PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%% Experimental Data: Responses to ethyl acetate
% What do we observe in ORN responses to odors? In the following figure, we plot the responses of ab3A to increasing pulses ethyl acetate, which are either presented in isolation (black), or on top of background (coloured plots). Each panel shows the experimentally determined input-output relationship for the neuron (dots), together with the best-fit Hill function. 

%%
% In the first panel, the exponent of the Hill Function is kept constant, while the other parameters are allowed to vary to fit the data. In the second panel, the offset of the Hill function is kept constant, and in the third, all parameters are allowed to vary. The legend shows the r-square between the data and the Hill function fit.

load('/local-data/DA-paper/dose-response/carlotta/dr_data.mat')
Lh= [];
L = {};
figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
for i = 1:3
	p(i) =  subplot(1,3,i); hold on
	% no background
	baseline = [mean(data(7).PID(2:70,:)) mean(data(1).PID(2:70,:))];
	x=[max(data(7).PID) max(data(1).PID)];
	y=[max(data(7).ORN) max(data(1).ORN)];
	scatter(p(i),x,y,'k')
	fo=optimset('MaxFunEvals',1000,'Display','none');
	hf0 = lsqcurvefit(@hill4,[max(y) 2 2 0],x(:),y(:),[max(y)/2 0 1 0],[2*max(y) max(y) 10 10],fo);
	Lh(1) = plot(p(i),sort(x),hill4(hf0,sort(x)),'k');
	L{1} = oval(rsquare(y,hill4(hf,x)),2);
end

backgrounds = [5 6];
title(p(1),'Divisive Only')
for i = 1:length(backgrounds)
	baseline = mean(data(backgrounds(i)).PID(2:70,:));
	x=max(data(backgrounds(i)).PID);
	y=max(data(backgrounds(i)).ORN);
	scatter(p(1),x,y,32,C(i,:))
	fo=optimset('MaxFunEvals',1000,'Display','none');
	hf = lsqcurvefit(@hill4,[max(y) hf0(2) 2 0],x(:),y(:),[max(y)/2 hf0(2) 0 0],[2*max(y) hf0(2)+1e-6 10 10],fo);
	Lh(i+1) = plot(p(1),sort(x),hill4(hf,sort(x)),'Color',C(i,:));
	div_only_hf = hf;
	L{i+1} = oval(rsquare(y,hill4(hf,x)),2);
end
xlabel(p(1),'Stimulus (V)')
ylabel(p(1),'Response (Hz)')
set(p(1),'XScale','log')
legend(p(1),Lh,L,'location','northwest')


title(p(2),'Subtractive Only')
for i = 1:length(backgrounds)
	baseline = mean(data(backgrounds(i)).PID(2:70,:));
	x=max(data(backgrounds(i)).PID);
	y=max(data(backgrounds(i)).ORN);
	scatter(p(2),x,y,32,C(i,:))
	fo=optimset('MaxFunEvals',1000,'Display','none');
	hf = lsqcurvefit(@hill4,[max(y) hf0(2) hf0(3) 0],x(:),y(:),[max(y)/2 0 hf0(3) 0],[2*max(y) 10 hf0(3)+1e-6 10],fo);
	Lh(i+1) = plot(p(2),sort(x),hill4(hf,sort(x)),'Color',C(i,:));
	sub_only_hf = hf;
	L{i+1} = oval(rsquare(y,hill4(hf,x)),2);
end
xlabel(p(2),'Stimulus (V)')
set(p(2),'XScale','log')
legend(p(2),Lh,L,'location','northwest')


title(p(3),'Mixed Model')
for i = 1:length(backgrounds)
	baseline = mean(data(backgrounds(i)).PID(2:70,:));
	x=max(data(backgrounds(i)).PID);
	y=max(data(backgrounds(i)).ORN);
	scatter(p(3),x,y,32,C(i,:))
	fo=optimset('MaxFunEvals',1000,'Display','none');
	hf = lsqcurvefit(@hill4,[max(y) 2 2 0],x(:),y(:),[max(y)/2 0 0 0],[2*max(y) max(y) 10 10],fo);
	Lh(i+1)=plot(p(3),sort(x),hill4(hf,sort(x)),'Color',C(i,:));
	mixed_hf = hf;
	L{i+1} = oval(rsquare(y,hill4(hf,x)),2);
end
xlabel(p(3),'Stimulus (V)')
set(p(3),'XScale','log')
legend(p(3),Lh,L,'location','northwest')


PrettyFig;
if being_published
	snapnow
	delete(gcf)
end




%%
% In the following figure, we repeat the analysis, but now plot responses as a function of stimulus - background. 

Lh= [];
L = {};
figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
for i = 1:3
	p(i) =  subplot(1,3,i); hold on
	% no background
	baseline = [mean(data(7).PID(2:70,:)) mean(data(1).PID(2:70,:))];
	x=[max(data(7).PID) max(data(1).PID)] - baseline;
	y=[max(data(7).ORN) max(data(1).ORN)];
	scatter(p(i),x,y,'k')
	fo=optimset('MaxFunEvals',1000,'Display','none');
	hf0 = lsqcurvefit(@hill4,[max(y) 2 2 0],x(:),y(:),[max(y)/2 0 1 0],[2*max(y) max(y) 10 10],fo);
	Lh(1) = plot(p(i),sort(x),hill4(hf0,sort(x)),'k');
	L{1} = oval(rsquare(y,hill4(hf,x)),2);
end

backgrounds = [5 6];
title(p(1),'Divisive Only')
for i = 1:length(backgrounds)
	baseline = mean(data(backgrounds(i)).PID(2:70,:));
	x=max(data(backgrounds(i)).PID)-baseline;
	y=max(data(backgrounds(i)).ORN);
	scatter(p(1),x,y,32,C(i,:))
	fo=optimset('MaxFunEvals',1000,'Display','none');
	hf = lsqcurvefit(@hill4,[max(y) hf0(2) 2 0],x(:),y(:),[max(y)/2 hf0(2) 0 0],[2*max(y) hf0(2)+1e-6 10 10],fo);
	Lh(i+1) = plot(p(1),sort(x),hill4(hf,sort(x)),'Color',C(i,:));
	div_only_hf = hf;
	L{i+1} = oval(rsquare(y,hill4(hf,x)),2);
end
xlabel(p(1),'Stimulus - Background (V)')
ylabel(p(1),'Response (Hz)')
set(p(1),'XScale','log')
legend(p(1),Lh,L,'location','northwest')


title(p(2),'Subtractive Only')
for i = 1:length(backgrounds)
	baseline = mean(data(backgrounds(i)).PID(2:70,:));
	x=max(data(backgrounds(i)).PID)-baseline;
	y=max(data(backgrounds(i)).ORN);
	scatter(p(2),x,y,32,C(i,:))
	fo=optimset('MaxFunEvals',1000,'Display','none');
	hf = lsqcurvefit(@hill4,[max(y) hf0(2) hf0(3) 0],x(:),y(:),[max(y)/2 0 hf0(3) 0],[2*max(y) 10 hf0(3)+1e-6 10],fo);
	Lh(i+1) = plot(p(2),sort(x),hill4(hf,sort(x)),'Color',C(i,:));
	sub_only_hf = hf;
	L{i+1} = oval(rsquare(y,hill4(hf,x)),2);
end
xlabel(p(2),'Stimulus - Background (V)')
set(p(2),'XScale','log')
legend(p(2),Lh,L,'location','northwest')


title(p(3),'Mixed Model')
for i = 1:length(backgrounds)
	baseline = mean(data(backgrounds(i)).PID(2:70,:));
	x=max(data(backgrounds(i)).PID)-baseline;
	y=max(data(backgrounds(i)).ORN);
	scatter(p(3),x,y,32,C(i,:))
	fo=optimset('MaxFunEvals',1000,'Display','none');
	hf = lsqcurvefit(@hill4,[max(y) 2 2 0],x(:),y(:),[max(y)/2 0 0 0],[2*max(y) max(y) 10 10],fo);
	Lh(i+1)=plot(p(3),sort(x),hill4(hf,sort(x)),'Color',C(i,:));
	mixed_hf = hf;
	L{i+1} = oval(rsquare(y,hill4(hf,x)),2);
end
xlabel(p(3),'Stimulus - Background (V)')
set(p(3),'XScale','log')
legend(p(3),Lh,L,'location','northwest')


PrettyFig;
if being_published
	snapnow
	delete(gcf)
end



% cache data
being_published = 0;
save(strcat(mfilename,'.mat'))

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



