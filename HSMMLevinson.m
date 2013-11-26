%%
clear all
close all
clc
%%
format long e
Y = load('longstan200.txt');
Y(:,1) = Y(:,1) - Y(1,1);
Y = sortrows(Y);
Y1 = Y(:,1);
Y2 = Y(:,2);
Y2 = medfilt1(Y2, 15);
L = length(Y2);
Y2 = Y2-min(Y2)+0.001;
%%
counter = 0;
epsilon = 1;
N = 4;%number of states
y = Y2(1:50000);%extract learning sequence
T = length(y);
ymax = max(y);

A = ones(N, N)./(N-1); %in HSMM, there's no transition to the same state
for i = 1:N
    A(i,i) = 0;
end

% A = ones(N,N)/N;

G = (1:N)*ymax/(N+1);
V = ones(N,1);

alpha = zeros(T, N);
beta = zeros(T,N);


mu = zeros(N,1); %this was called d in the previous codes and the paper, but now we need d for state duration
sigma = zeros(N,1);

Ahat = zeros(N,N);
muhat = zeros(N,1); %cf mu
sigmahat = zeros(N,1);


step = 0.001; %0.001 seems to be the resolution of the measured delay anyway
x = 0:step:(floor(ymax/step)+1)*step; 
b = zeros(length(x),N);

c = ones(T,1);

%for HSMM:
param = zeros(N,2);
param(1,:) = [4.64 1.328]; %how to choose those ?
param(2,:) = [3.08 0.831]; %this is from the preparation....
param(3,:) = [2 0.5]; %au pif....
D = 4000; %maximum number of timesteps we can stay in one state %4000 roughly...
%1000 might be too much computation intensive
Ddist = zeros(N,D); %this is the p of the paper

% %this is for a HMM
% D = 1;
% Ddist(:,1)=1;

xd = 1:1:D ;
for i = 1:N
    pd = makedist('logn', 'mu',param(i,1),'sigma',param(i,2) );
    %pd = makedist('Poisson','lambda', lambda);
    Ddist(i,:) = (pdf(pd,xd));
end

%those qre the parameters of the duration distribution. Needs to be
%initialized properly.
shape = zeros(N,1);
scale = zeros(N,1);
shapehat = shape;
scalehat = scale;



alpha(D,1) = 1;
beta(T-D+1,:) = ones(1,N); 
alpha0 = alpha;
beta0 = beta;
%%
%while(counter <30)
tic;
    counter= counter+1;
    display(counter)
    display('Generate distributions');
    for i=1:N
        %observations:
        pd = makedist('Gamma','a',G(i),'b',V(i));
        b(:,i) = (pdf(pd,x))';
    end
    
    %
    display('Compute alphas');
    alpha = alpha0;
    for t=D+1:T-D
        display(t);
        for j=1:N
            for d=1:D 
                for i=1:N
                    term = alpha(t-d,i)*A(i,j)*Ddist(j,d)*prod(c(t-d+1:t-1));
                    term =term*prod(b(round(y(t-d+1:t)/step)+1,j));
                    alpha(t,j) = alpha(t,j)+term;
                end
            end
        end
        c(t) = 1./sum(alpha(t,:));
        alpha(t,:) = alpha(t,:).*c(t);
    end
    
    
    display('Compute betas');
    beta = beta0;
    for t=T-D:-1:D+1
        display(t);
        for i=1:N
            for d=1:D
                for j=1:N
                    term = A(i,j)*Ddist(j,d)*prod(c(t+1:t+d-1))*beta(t+d,j);
                    term = term*prod(b(round(y(t+1:t+d)/step)+1,j));
                    beta(t,i) = beta(t,i)+term;
                end
            end
        end
        beta(t,:) = beta(t,:).*c(t);
    end



    display('Compute A');
    for i=1:N
        for j=1:N
            display(i);
            display(j);
            num = 0;
            for t=D+1:T-D
                for d=1:D
                    term = A(i,j)*Ddist(j,d)*beta(t+d,j)*prod(c(t+1:t+d-1));
                    term = term*alpha(t,i)*prod(b(round(y(t+1:t+d)/step)+1,j));
                    num = num + term;
                end
            end
            den = sum(alpha(D+1:T-D,i).*beta(D+1:T-D,i)./c(D+1:T-D));  
            Ahat(i,j) = num/den;
        end
    end








    display('Compute mu and sigma');
    for j=1:N
        display(j);
        nummu = 0;
        numsig = 0;
        den = 0;
        for t = D+1:T-D
            for d = 1:D
                for i=1:N
                    term = alpha(t-d,i)*A(i,j)*Ddist(j,d)*prod(c(t-d+1:t-1));
                    term = term*prod(b(round(y(t-d+1:t)/step)+1,j));
                    nummu = nummu+term*sum(y(t-d+1:t));
                    numsig = numsig+term*sum((y(t-d+1:t)-mu(j)).^2);
                    den = den+term; %den+term or den+term*d ???
                end
            end
        end
        muhat(j) = nummu/den;
        sigmahat(j) = numsig/den;
    end

    
    display('Compute duration parameters');
    %Update duration parameters:
    %scale:
    for j=1:N
        num = 0;
        for t=D+1:T-D
            display(j);
            display(t);
            for d=1:D
                for i=1:N
                    term = d*alpha(t-d,i)*A(i,j)*Ddist(j,d)*prod(c(t-d+1:t-1));
                    term = term*beta(t,j)*prod(b(round(y(t-d+1:t)/step)+1,j));
                    num = num + term;
                end
            end
        end
        den = shape(j)*sum(alpha(D+1:T-D,j).*beta(D+1:T-D,j)./c(D+1:T-D));
        scalehat(j) = num/den;
    end
    
    %shape:
    
    mu = muhat;
    sigma = sigmahat;
    A = Ahat;
    %convert to gamma dist parameters:
    V = sigma./mu;
    G = mu.^2./sigma;
    
    scale = scalehat;
    shape = shapehat;
    rate = 1./scale;
    
    L = -sum(log(c));
    display(L);
    
    
    
    display('Compute duration distribution parameters');
    C = zeros(N,1);
    digamma = psi(shape);
    polygamma=psi(1,shape);
    for j=1:N
        num = 0;
        for t=D+1:T-D
            for d=1:D
                for i=1:N
                    term = log(rate(j)*d)*alpha(t-d,i)*A(i,j)*Ddist(j,d);
                    term = term*prod(b(round(y(t-d+1:t)/step)+1,j));
                    term = term*beta(t,j)*prod(c(t-d+1:t-1));
                    num = num + term;
                end
            end
        end
        den = sum(alpha(D+1:T-D,j).*beta(D+1:T-D, j)./c(D+1:T-D));
        C(j) = num/den; 
        shapehat(j) = shapehat(j) -(digamma(j)-C(j))/(polygamma(j));
    end
    
    
    
    
    
    display('done.')
    
    
    
    
    
%end






















