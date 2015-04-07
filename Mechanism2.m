% Mechanism2.m
% in this document, we attempt to understand the mechanism behind gain adaptation
% 
% created by Srinivas Gorur-Shandilya at 1:15 , 06 April 2015. Contact me at http://srinivas.gs/contact/
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


%% What is the mechanism behind gain adaptation?
% Two (non-exclusive) possibilities are that the gain is controlled at the receptor level, and that the gain is controlled at the firing machinery. To disambiguate the two, we record the neuron's responses to a flickering odour stimulus. We then repeat the experiment, but also drive the neuron optically through light activating ReaChR channels. 

%% LEDs can deliver lots of power through the objective
% We deliver light through the objective, as shown below:
%
% <</Users/sigbhu/code/da/images/led.jpg>>
%

%%
% The following figure shows the light levels measured at the fly position vs the voltage driving the LED circuit. 

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
load('LED_calib.mat')
plot(LED.x,LED.y,'k+')
xlabel('Driving Voltage (V)')
ylabel('Power @ 627nm (mW)')
clear p
p.A= 2.7288;
p.k= 1.9863;
p.n= 2.8827;
l=plot(LED.x,hill(LED.x,p),'r')
legend(l,'Hill Fit','location','southeast')

PrettyFig();

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
