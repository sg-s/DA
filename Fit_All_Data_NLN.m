% fit all of carlotta's data with the NLN model
% 
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
load('/local-data/DA-paper/data.mat')

if ~exist('NLNParam.mat','file')

	NLNParam = zeros(4,21);

	for i = 2:21
		disp(data(i).original_name);
		clear d
		d.stimulus = data(i).PID;
		d.stimulus(d.stimulus<0) = 0;
		d.response = data(i).ORN;
		[~,x] = FitNLNModel(d,[],.99);
		NLNParam(:,i) = x;
	end

	save('NLNParam.mat','NLNParam')

else
	load('NLNParam.mat','NLNParam')

end

% check that all the fits are good. 

for i = 2:21
	figure, hold on

	% solve the NLN model
	data(i).PID(data(i).PID<0) =0 ;
	[Rguess] = SolveNLNModel2(NLNParam(:,i),data(i).PID,data(i).ORN);

	plot(data(i).time,data(i).ORN,'k')
	plot(data(i).time,Rguess,'r')

	title(data(i).original_name)

	PrettyFig;

end