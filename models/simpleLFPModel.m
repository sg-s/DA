% simpleLFPModel.m
% inspired by the DA model, but attempts to get the LFP timescales right
% 
% created by Srinivas Gorur-Shandilya at 11:19 , 26 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [R] = simpleLFPModel(S,p)

% parameters
p.tau; % timescale
p.m; % exponent 

p.A;
p.B;

lb.tau = eps;
lb.m = .1; 


if size(S,1) < size(S,2)
	S = S';
end

R = 0*S;

for j = 1:size(S,2)
	S(:,j) = filter(ones(10,1),10,S(:,j));
	for i = 2:length(S)
		fx1 = p.A*(S(i-1,j).^(-1 - p.m));
		fx2 = p.B*(S(i-1,j).^(-p.m))*R(i-1,j);
		fx = (fx1-fx2)/p.tau;

		R(i,j) = fx + R(i-1,j);
	end

	% normalise
	R(:,j) = R(:,j) - mean(R(3e3:5e3,j));
	R(:,j) = -R(:,j);
	R(:,j) = R(:,j)/max(R(5e3:6e3,j));
end

