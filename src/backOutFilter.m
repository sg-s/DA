% backOutFilter.m
% Baccus-Meister-style filter extraction using FFT
% modified from Damon's original code
%  
% 


function [K, varargout]= backOutFilter(x,R,varargin)

% warning: x and R are row vectors, not column vectors
x=x(:)' - mean(x); R=R(:)' - mean(R);


% defaults
filter_length = 700;
offset = 100;
filter_buffer = 100;
EPS = 0;
ZEROPAD = true;
DOXCORR = false;

if ~nargin
    help functionName
    return
else
    if iseven(length(varargin))
        for ii = 1:2:length(varargin)-1
            temp = varargin{ii};
            if ischar(temp)
                eval(strcat(temp,'=varargin{ii+1};'));
            end
        end
    else
        error('Inputs need to be name value pairs')
    end
end

% handle the offset
if offset ~= 0 
    x = [x NaN(1,offset+round(filter_buffer/2))]; 
    R = [NaN(1,offset+round(filter_buffer/2)) R];
end
rm_this = isnan(x) | isnan(R);
x(rm_this) = []; R(rm_this) = [];
N = filter_length + filter_buffer;

interval = 5; % so you aren't in phase with 60 hz pick up
spacing=1:interval:length(R)-N;

if ~ZEROPAD
    xw=zeros(length(spacing),N);
    Rw=xw;
else
    xw=zeros(length(spacing),2*N);
    Rw=xw;
end


for i=1:length(spacing)
    if ~ZEROPAD
        xw(i,:)=fft(x(spacing(i):spacing(i)+N-1));
        Rw(i,:)=fft(R(spacing(i):spacing(i)+N-1));
    else
        xw(i,:)=fft([x(spacing(i):spacing(i)+N-1),zeros(1,N)]);
        Rw(i,:)=fft([R(spacing(i):spacing(i)+N-1),zeros(1,N)]);
    end
end


num=mean((Rw).*conj(xw),1);
ssdenom = sum(abs(mean(xw.*conj(xw),1)).^2);
denom=mean(xw.*conj(xw),1)+EPS*abs(mean(mean(xw.*conj(xw),1)));
ssdenom2 = sum(abs(denom).^2);
denom = denom * (ssdenom/ssdenom2)^0.5; % equalize power with EPS!

if ~DOXCORR
    Kw = num./denom; % DAC 081216 -- conjugate stimulus, as in baccus linear filter
else
    Kw = num/std(x); % don't put in the stimulus stats
end

K=ifft(Kw);
if ZEROPAD
    K=K(1:end/2); %% IS THIS THE CORRECT HALF? think so.
end

Rt=filter(K,1,x);
C=corr([R(:),Rt(:)]);
varargout{1}=C(1,2);
varargout{2}=Rt;

% chop off the buffer
K = K(round(filter_buffer/2)+1:end-round(filter_buffer/2));


% normlaise filter
K = K/sqrt(sum(K.^2));

