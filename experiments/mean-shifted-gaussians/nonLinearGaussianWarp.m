% nonLinearGaussianWarp
% this function is meant to be run by fitModel2Data
% its input is a control signal that goes to the MFC
% its output is a 3-element long vector: [mean(predictedPID) std(predictedPID) 100*min(predictedPID) 100*max(predictedPID) 1-r2(Gaussianfit)]
% it assumes that there is are globals called K and ff (together, the LN model)
% this is how it works:
% it accepts the stimulus, and warps it using a nonlinear function defined by the parameters structure
% after warping, it uses the LN model defined in the global workspace to make a prediction. 
% and then calculates the output as defined above.
%
% created by Srinivas Gorur-Shandilya at 1:55 , 10 September 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function r = nonLinearGaussianWarp(stim,p)

global K
global ff

% define parameters
p.c1;
p.c2;
p.c3;
p.c4;
p.c5;
p.c6;

% warp the control signal
warped_stim = p.c1*stim.^5 + p.c2*stim.^4 + p.c3*stim.^3 + p.c4*stim.^2 + p.c5*stim + p.c6; 

% prevent ridiculous values
mfc_min = 5/200;
mfc_max = 5;
warped_stim(warped_stim < mfc_min) = mfc_min;
warped_stim(warped_stim > mfc_max) = mfc_max;


% make a prediction of the PID
PID_pred = filter(K,1,warped_stim);
PID_pred = ff(PID_pred);



r = hist(PID_pred,0:0.01:5);
r = r/max(r);
