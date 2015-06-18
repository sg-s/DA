function [K, varargout]= Back_out_1dfilter_new(x,R,varargin);
% code from Damon, found on the wiki...
% attempt to do the baccus-style filtering to get the two dimensional
% filter back out...

% 090429 -- this is not working all of a sudden. perhaps just not working
% with binary inputs? doesn't work well with width 2 delta functions...

x=x(:)' - mean(x); R=R(:)' - mean(R);

if ~length(varargin)
    N=100; % size of filter to look for...
else
    N=varargin{1};
end
DOXCORR = 1; % don't do by default
ZEROPAD = 1; % 1 is default
if length(varargin)>1
    if varargin{2}==-1
        DOXCORR=1;
		EPS = 0;
    else
        EPS = varargin{2};
    end
else
    EPS = 0;
end

%interval = 3; % so you aren't in phase with 60 hz pick up
interval = 1;
spacing=1:interval:length(R)-N;

if ZEROPAD
    xw=zeros(length(spacing),2*N);
    Rw=xw;
else
    xw=zeros(length(spacing),N);
    Rw=xw;
end


for i=1:length(spacing)
    if ZEROPAD
        xw(i,:)=fft([x(spacing(i):spacing(i)+N-1),zeros(1,N)]);
        Rw(i,:)=fft([R(spacing(i):spacing(i)+N-1),zeros(1,N)]);
        
    else
        xw(i,:)=fft(x(spacing(i):spacing(i)+N-1));
        Rw(i,:)=fft(R(spacing(i):spacing(i)+N-1));
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


