% nonLinearWarpMSG.m
% this function nonlinearly warps the control signals of a MFC so that PIDs that come out of it are more Gaussian distributed. 
% 
% the way it works is this:
% 0. first, you need to generate some data with real PID measurements using some ansatz controls
% then, for each control paradigms, nonLinearWarpMSG:
% 1. builds a LN model from the control signals to the PID 
% 2. picks a nonlinear function to be applied to the control signals so that you minimise:
% a) the difference in means from the warped PIDs and the original PIDs
% b) the difference in stds from " " " "
% c) 1 - r^2(a guassian fit to the histogram of the warped PIDs)
%
% created by Srinivas Gorur-Shandilya at 10:06 , 10 September 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. x
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

% assume data is loaded in workspace, duplicate ControlParadigms and data
warped_data = data;
warped_ControlParadigm = ControlParadigm;

% some timing and other parameters
a = 15e4;
z = 55e4;
mfc_output = 1;

global K
global ff
LN_model = struct;

for i = 1:length(data)
	if ~isempty(data(i).PID)
		x = ControlParadigm(i).Outputs(mfc_output,a:10:z);
		y = mean2(data(i).PID(:,a:10:z));

		disp('Extracting Linear filter...')
		% fit a LN model to this
		LN_model(i).K = fitFilter2Data(x,y,'filter_length',500,'reg',1);

		pred = filter(LN_model(i).K,1,x);
		disp('Fitting nonlinearity...')
		LN_model(i).ff = fit(pred(:),y(:),'poly5');

		disp('DONE. Model fit quality is:')
		disp(rsquare(LN_model(i).ff(pred), y));

		disp('Finding best warp of data...')
		clear d
		d.stimulus = x;
		d.response = normpdf(0:0.01:5,mean(y),std(y));
		d.response = d.response/max(d.response);

		K = LN_model(i).K;
		ff = LN_model(i).ff;

		p = fitModel2Data(@nonLinearGaussianWarp,d,'UseParallel',false);

		% now warp the control paradigms
		stim = ControlParadigm(i).Outputs(1,:);
		temp = p.c1*stim.^5 + p.c2*stim.^4 + p.c3*stim.^3 + p.c4*stim.^2 + p.c5*stim + p.c6; 
		warped_ControlParadigm(i).Outputs(mfc_output,:) = temp;

	end

end

ControlParadigm = warped_ControlParadigm;
save('warped_MSG_Kontroller_paradigm.mat','ControlParadigm')




