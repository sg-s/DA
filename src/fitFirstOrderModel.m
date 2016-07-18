%% fitFirstOrderModel.m
% fits a first order model:
% 0 = K0 x s(t) + (1 + K1 x s(t)) r(t)
% where x is a convolution
% to data specified by s(t) and r(t)
% and K0 and K1 are filters that are fit non-parametrically 


function [K0,K1] = fitFirstOrderModel(S,R,varargin)



% options and defaults
options.filter_length = 500;
options.reg = 0;
options.whiten = 'true';


if nargout && ~nargin 
	varargout{1} = options;
    return
end

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
elseif isstruct(varargin{1})
	% should be OK...
	options = varargin{1};
else
	error('Inputs need to be name value pairs')
end


assert(isvector(S),'First argument must be a vector')
assert(isvector(R),'2nd argument should be a vector')
assert(length(S) == length(R),'first two arguments should be of equal length')

S = S(:);
R = R(:);

only_these_points = options.filter_length+1:length(R);

R = R(only_these_points);

% chop up the stimulus into blocks  
stim = zeros(options.filter_length,length(only_these_points));

for i = 1:length(only_these_points)
	stim(:,i) = S(only_these_points(i):-1:only_these_points(i)-options.filter_length+1);
end

% also multiply each row by the response and save that too
stim_r = stim;
for i = 1:size(stim_r,1)
	stim_r(i,:) = stim_r(i,:).*R';
end

% combine these two
shat = [stim; stim_r];

% extract filters, based on regularisation
if ~options.whiten 
	% disp('No whitening. Directly computing filters...')
	K = R'/shat;

else
	% disp('Whitening stimulus...')

	% whiten and regularize stimulus separately for each filter
	C1 = stim*stim';
	MeanEigenValue = trace(C1)/length(C1); % cheat; this is the same as mean(eig(C1))
	r = options.reg*MeanEigenValue;
	C1 = (C1 + r*eye(length(C1)))*trace(C1)/(trace(C1) + options.reg*length(C1));

	C2 = stim_r*stim_r';
	MeanEigenValue = trace(C2)/length(C2); % cheat; this is the same as mean(eig(C2))
	r = options.reg*MeanEigenValue;
	C2 = (C2 + r*eye(length(C2)))*trace(C2)/(trace(C2) + options.reg*length(C2));

	C = shat*shat';
	MeanEigenValue = trace(C)/length(C); % cheat; this is the same as mean(eig(C))
	r = options.reg*MeanEigenValue;
	C = (C + r*eye(length(C)))*trace(C)/(trace(C) + options.reg*length(C));
	C(1:length(C1),1:length(C1)) = C1;
	C(length(C1)+1:end,length(C1)+1:end) = C2;
    

    K = C\(shat*R);

end

K0 = K(1:length(K)/2);
K1 = K(length(K)/2+1:end);




