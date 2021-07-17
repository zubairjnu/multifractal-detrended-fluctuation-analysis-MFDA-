% actual code name GFGU_MFDMA_1D.m
%
function [n,Fq,tau,alpha,f] = F_ALPHA(x,n_min,n_max,N,theta,q)
%
% The procedure [n,Fq,tau,alpha,f]=GFGU_MFDMA_1D(x,n_min,n_max,N,theta,q) is
% used to calculate the multifractal properties of one-dimensional time series.
%
% Input:
%   x: the time series we considered
%   n_min: the lower bound of the segment size n
%   n_max: the upper bound of the segment size n
%   N: the length of n, that is, the data points in the plot of Fq VS n
%   theta: the position parameter of the moving window
%   q: multifractal order

if size(x,2) == 1
    x = x';
end

M = length(x);
MIN = log10(n_min);
MAX = log10(n_max);
n = (unique(round(logspace(MIN,MAX,N))))';

% Construct the cumulative sum y
y = cumsum(x);

for i = 1:length(n)
    lgth = n(i,1);

    % Calculate the moving average function \widetilde{y}
    y1 = zeros(1,M-lgth+1);
    for j = 1:M-lgth+1
        y1(j) = mean(y(j:j+lgth-1));
    end

    % Determine the residual e
    e=y(max(1,floor(lgth*(1-theta))):max(1,floor(lgth*(1-theta)))+length(y1)-1)-y1;

    % Estimate the root-mean-square function F
    for k=1:floor(length(e)/lgth)
        F{i}(k)=sqrt(mean(e((k-1)*lgth+1:k*lgth).^2));
    end
end


% Calculate the q-th order overall fluctuation function Fq
for i=1:length(q)
    for j=1:length(F)
        f=F{j};
        if q(i) == 0
            Fq(j,i)=exp(0.5*mean(log(f.^2)));
        else
            Fq(j,i)=(mean(f.^q(i)))^(1/q(i));
        end
    end
end


% Calculate the multifractal scaling exponent tau(q)
for i=1:size(Fq,2)
    fq=Fq(:,i);
    r=regstats(log(fq),log(n),'linear',{'tstat'});
    k=r.tstat.beta(2);
    h(i,1)=k;
end
tau=h.*q'-1;


% Calculate the singularity strength function alpha(q) and spectrum f(alpha)
dx=7;
dx=fix((dx-1)/2);
for i=dx+1:length(tau)-dx
    xx=q(i-dx:i+dx);
    yy=tau(i-dx:i+dx);
    r=regstats(yy,xx,'linear',{'tstat'});
    alpha(i,1)=r.tstat.beta(2);
end
alpha=alpha(dx+1:end);
f=q(dx+1:end-dx)'.*alpha-tau(dx+1:end-dx);
