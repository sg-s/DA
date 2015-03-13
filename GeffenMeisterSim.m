% GeffenMeisterSim
% 
% created by Srinivas Gorur-Shandilya at 12:53 , 13 March 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [] = GeffenMeisterSim(p)

% check if the figure has been created
DASim_figure = [];
try GMSim_figure = evalin('base', 'GMSim_figure');
	GMSim_plot = evalin('base', 'GMSim_plot');
catch
	GMSim_figure = figure('Units','pixels','Position',[10 100 750 700],'Name','Geffen-Meister-Sim','IntegerHandle','off','CloseRequestFcn',@closess,'Color','w'); hold on
	clf;
	GMSim_plot(1) = axes('Units','pixels','Position',[50 450 650 200]);
	GMSim_plot(2) = axes('Units','pixels','Position',[50 100 300 300]);
	GMSim_plot(3) = axes('Units','pixels','Position',[450 100 200 200]);

	% save to the base workspace
	assignin('base','GMSim_figure',GMSim_figure)
	assignin('base','GMSim_plot',GMSim_plot)
end

% some hard coded parameters
T = 60;
dt = 1e-2;
time = dt*(1:floor(T/dt));
time = time(:);
marker_size=10;
marker_size2=24;
font_size=20;

% bounds

lb.ntrials = 2;
ub.ntrials = 20;
lb.w = 0;
ub.w = 2;
lb.what = 0;
ub.what = 2;
ub.noise = 3;
lb.A = 0;
ub.A = 2;

if p.ntrials  < 2
	p.ntrials =2;
end

data = NaN(floor(T/dt),round(p.ntrials));
for i = 1:round(p.ntrials)
	data(:,i) = sin(time*p.w) + p.noise*randn(length(time),1);
end
plot(GMSim_plot(1),time,data,'k')

% plot the model predictions
pred = sin(time*p.what)*p.A;
hold(GMSim_plot(1),'on')
plot(GMSim_plot(1),time,pred,'r')
r2 = rsquare(pred,mean2(data));
title(GMSim_plot(1),strcat('r^2=',oval(r2)))
hold(GMSim_plot(1),'off')


[qx, qy] = GeffenMeister(data,pred);
plot(GMSim_plot(2),qx,qy,'r+')
hold(GMSim_plot(2),'on')
m = max([qx qy]*1.5);
set(GMSim_plot(2),'XLim',[0 m],'YLim',[0 m])
plot(GMSim_plot(2),[0 m],[0 m],'k--')
xlabel(GMSim_plot(2),'(P_{S}/P_{N})^{1/2}','interpreter','tex')
ylabel(GMSim_plot(2),'(P_{S}/P_{R})^{1/2}','interpreter','tex')
hold(GMSim_plot(2),'off')

plot(GMSim_plot(3),p.w,1,'ko')
hold(GMSim_plot(3),'on')
plot(GMSim_plot(3),p.what,p.A,'r+')
xlabel(GMSim_plot(3),'Frequency')
ylabel(GMSim_plot(3),'Amplitude')
set(GMSim_plot(3),'XLim',[0 2],'YLim',[0 2])
hold(GMSim_plot(3),'off')

function [] = closess(~,~)
	delete(GMSim_figure)
	evalin('base','clear GMSim_figure')
	evalin('base','clear GMSim_plot')

end

end