%% Plot Interest Rates Against Time
clear all

%load US10Y.mat
load yielddata.mat


%% Find Startdate and enddate
[x,y] = find(tGermany==Germanystartenddate9m(2,1));
% x=1099;
% y=1;
startdate=tGermany(x,y);
[a,y] = find(tGermany==Germanystartenddate9m(2,2));             
%a=3243;
y=1;
enddate=tGermany(a,y);

t=tGermany(x:a);
data=Germany(x:a);
%% Plot yields against time

% t=t+693960;
t1 = t+270;
%data1=UK(x+270:a+270);
figure
plot(t1,data);
datetick('x','yyyy-mm')
title('10Y Treasury Yield');
xlabel('Date');
ylabel('Yield (%)');

T=datestr(t);
%% 

logR=log(data);
%logR=US;
YieldTimes = yearfrac(t(1), t);

seasonMatrix = @(t) [sin(2.*pi.*t) cos(2.*pi.*t) sin(4.*pi.*t) cos(4.*pi.*t) t ones(size(t, 1), 1)];
%seasonMatrix = @(t) [ t ones(size(t, 1), 1)];
C = seasonMatrix(YieldTimes);
seasonParam = C\logR;

%%
figure
subplot(2, 1, 1);
plot(t, logR);
datetick('x','yyyy-mm');
title('log(r) and deterministic fluctuations');
xlabel('Date');
ylabel('log(r)');
hold on;
plot(t, C*seasonParam, 'r');
hold off;
legend('log(r)', 'deterministic fluctuations');
X = logR-C*seasonParam;
subplot(2, 1, 2);
plot(t, X);
datetick('x','yyyy-mm');
title('log(R) with Deterministic Fluctuations Removed');
xlabel('Date');
ylabel('log(R)');

%%
% Prices at t, X(t)
Pt = X(2:end);
%Pt = logR(2:end);

% Prices at t-1, X(t-1)
 Pt_1 = X(1:end-1);
%Pt_1 = logR(1:end-1);

% Discretization for daily prices
dt = 1/270;

% PDF for discretized model
mrjpdf = @(Pt, a, phi, mu_J, sigmaSq, sigmaSq_J, lambda) ...
    lambda.*exp((-(Pt-a-phi.*Pt_1-mu_J).^2)./ (2.*(sigmaSq+sigmaSq_J))).* (1/sqrt(2.*pi.*(sigmaSq+sigmaSq_J))) + (1-lambda).*exp((-(Pt-a-phi.*Pt_1).^2)/(2.*sigmaSq)).* (1/sqrt(2.*pi.*sigmaSq));

%mrjpdf = @(Pt, phi, mu_J, sigmaSq, sigmaSq_J, lambda) ...
   % lambda.*exp((-(Pt-phi.*Pt_1-mu_J).^2)./ (2.*(sigmaSq+sigmaSq_J))).* (1/sqrt(2.*pi.*(sigmaSq+sigmaSq_J))) + (1-lambda).*exp((-(Pt-phi.*Pt_1).^2)/(2.*sigmaSq)).* (1/sqrt(2.*pi.*sigmaSq));


% Constraints: 
% phi < 1 (k > 0)
% sigmaSq > 0
% sigmaSq_J > 0
% 0 <= lambda <= 1
lb = [-Inf -Inf -Inf 0 0 0];
ub = [Inf 1 Inf Inf Inf 1];

% Initial values
x0 = [0 0 0 var(X) var(X) 0.5];

% Solve maximum likelihood
params = mle(Pt,'pdf',mrjpdf,'start',x0,'lowerbound',lb,'upperbound',ub,'optimfun','fmincon');

% Obtain calibrated parameters
alpha = params(1)/dt
 kappa = params(2)/dt
 mu_J = params(3)
 sigma = sqrt(params(4)/dt)
 sigma_J = sqrt(params(5))
 lambda = params(6)/dt

% kappa = params(1)/dt
% mu_J = params(2)
% sigma = sqrt(params(3)/dt)
% sigma_J = sqrt(params(4))
% lambda = params(5)/dt


%%
rng default;

% Simulate for about 2 years
nPeriods = 365+20;
nTrials = 10000;
n1 = randn(nPeriods,nTrials);
n2 = randn(nPeriods, nTrials);
j = binornd(1, lambda*dt, nPeriods, nTrials);
SimPrices = zeros(nPeriods, nTrials);
SimPrices(1,:) = X(end);
for i=2:nPeriods
    SimPrices(i,:) = alpha*dt + (1-kappa*dt)*SimPrices(i-1,:) + ...
                sigma*sqrt(dt)*n1(i,:) + j(i,:).*(mu_J+sigma_J*n2(i,:));
%                SimPrices(i,:) = (1-kappa*dt)*SimPrices(i-1,:) + ...
%                  sigma*sqrt(dt)*n1(i,:) + j(i,:).*(mu_J+sigma_J*n2(i,:));
end

% Add back seasonality
SimPriceDates = daysadd(t(end),0:nPeriods-1);
SimPriceTimes = yearfrac(t(1), SimPriceDates);
CSim = seasonMatrix(SimPriceTimes);
logSimPrices = SimPrices + repmat(CSim*seasonParam,1,nTrials);
%logSimPrices = SimPrices;

% Plot logarithm of Prices and simulated logarithm of Prices
figure
subplot(2, 1, 1);
plot(t, logR);
hold on;
plot(SimPriceDates(2:end), logSimPrices(2:end,1), 'red');
datetick('x','yyyy-mm');
 seasonLine = seasonMatrix([t; SimPriceTimes(2:end)])*seasonParam;
% plot([t; SimPriceDates(2:end)], seasonLine, 'green');
 %datetick('x','yyyy');
hold off;
title('Actual log(R) and Simulated log(price)');
xlabel('Date');
ylabel('log(R)');
legend('market', 'simulation');

% Plot prices and simulated prices
PricesSim = exp(logSimPrices);
subplot(2, 1, 2);
plot(t, data);
hold on;
plot(SimPriceDates, PricesSim(:,1), 'red');
hold off;
datetick('x','yyyy-mm');
title('Actual Prices and Simulated Prices');
xlabel('Date');
ylabel('Yield(%)');
legend('market', 'simulation');