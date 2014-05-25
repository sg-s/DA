% V2f.m
% accepts a matrix that is a set of raw voltage trace from ORN recordings, and outputs a firing rate. 
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

% a function that converts a matrix of spiketimes to a firing rate vector,
% or matrix.
% created by Srinivas Gorur-Shandilya at 15:38 , 19 July 2011. Contact me
% at http://srinivas.gs/contact/
function [fA,fB] = V2f(V,dt,t_on,gw,sliding_step)

% validate inputs
if nargin < 1
	error('No inputs?');
end	

if nargin < 2
	dt = 1e-4;
end	

if nargin < 3
	t_on = 1;
end	

if nargin < 4
	gw = 0.01;
end

if nargin < 5
	sliding_step = 0.003;
end

time =  dt:dt:dt*length(V);
x = sliding_step:sliding_step:max(time);
fA = zeros(1,length(x));
fB = zeros(1,length(x));

if isvector(V)
	% only one trial
	% find spikes
	[A,B] = V2Spike(V,dt,t_on);

	% add up the Gaussians
	for j = 1:length(A)
	    fA = fA + (normpdf(x,A(j)*dt,gw));   
	end


	for j = 1:length(B)
	    fB = fB + (normpdf(x,B(j)*dt,gw)); 
	end


else
	% many trials. pool all data
	for i = 1:size(V,1)

		[A,B] = V2Spike(V(i,:),dt,t_on);

		% add up the Gaussians
		for j = 1:length(A)
		    fA = fA + (normpdf(x,A(j)*dt,gw));   
		end


		for j = 1:length(B)
		    fB = fB + (normpdf(x,B(j)*dt,gw)); 
		end
	end

	% divide by number of trials
	fA = fA/size(V,1);
	fB = fB/size(V,1);
end

% debug
figure, hold on
plot(fB)
plot(fA,'r')


