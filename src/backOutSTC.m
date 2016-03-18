% backOutSTC.m
% wrapper function around Pillow's iSTAC code
% 
% created by Srinivas Gorur-Shandilya at 2:32 , 09 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [K1,K2,STA,KLcontributed] = backOutSTC(S,R,varargin)

% defaults
options.filter_length = 200;
options.n_dims = 10;
options.eigvalthresh = .05; % eigenvalue cutoff threshold (for pruning dims from raw stimulus)

% validate and accept options
if iseven(length(varargin))
	for ii = 1:2:length(varargin)-1
	temp = varargin{ii};
    if ischar(temp)
    	if ~any(find(strcmp(temp,fieldnames(options))))
    		disp(['Unknown option: ' temp])
    		disp('The allowed options are:')
    		disp(fieldnames(options))
    		error('UNKNOWN OPTION')
    	else
    		options = setfield(options,temp,varargin{ii+1});
    	end
    end
end
else
	error('Inputs need to be name value pairs')
end

% back out filters
[STA,stc,rawmu,rawcov] = simpleSTC(S,R,options.filter_length);


[vecs, vals] = compiSTAC(STA, stc, rawmu, rawcov, options.n_dims,options.eigvalthresh);
KLcontributed = [vals(1); diff(vals)];

K1 = vecs(:,1);
K2 = vecs(:,2);