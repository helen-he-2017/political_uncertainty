clc;clear all;close all;

%% import data
load yielddata.mat
rdata = US';
%% Input Values
nsample = 10;
[x,y] = find(tUS==USstartenddate1(1,1));
[a,y] = find(tUS==USstartenddate1(1,2));
%x=2107;
x=x+15;
startdate=x; %input start date
enddate=a;
%sigma = var(rdata(startdate:enddate))^0.5 %input total volatility (1 number)
sigma=0.05; %input parameters
T =length(rdata(startdate:enddate)); %input dates
nsimu=4000; %input simulations
%%
sigmagsample = linspace(0.1,1,nsample); %create a set of volatility
%% Simulations


rt=zeros(1,T); %store simulation result
rsimu=zeros(nsample,T); %store path simulation results
rsimu(:,1)= rdata(startdate) ; %input initial value

gtbar=zeros(nsample,T); %store posterior means over time
% gtbar(:,1)=rdata(startdate);
gtbar(:,1)=0;
sigmatbar=zeros(nsample,T); %store posterior volatility over time
sigmatbar(:,1)=1;
%stdsimu=zeros(1,nsample);

theta=zeros(nsample,T);
theta(:,1)=randn(nsample,1);
for k=1:nsample
    n1 = randn(nsimu,T);
    n2 = randn(1, T);
    n3 = randn(1,T);
    for t = 2:T
%         A=gtbar(:,t-1)+n1*sigmatbar(:,t-1);
        A=theta(k,t-1)+n1;
%         C=gtbar(:,t-1)+cumsum(n1)*sigmatbar(:,t-1);
        gtbar(k,t)=1/var(A(:,t-1))/(1/var(A(:,t-1))+1/sigmatbar(k,t-1))*1/nsimu*A(end,t);
        sigmatbar(k,t)=1/(1/sigmatbar(k,t-1)^2+1/var(A(:,t-1)));
        theta(k,t)= gtbar(k,t)+ sigmatbar(k,t).* n3(:,t);
        rsimu(k,t)=rsimu(k,t-1)+sigma*n2(:,t)+theta(k,t);
        %rsimu(:,t)=rsimu(:,t-1)+sigma*n1(:,t)+sigmatbar(:,t).*n2(:,t);
    end
end

    B=cumsum(rsimu);
    rt(1,:)=1/nsample*B(end,:); 
   rresult=cumsum((rt(1,:)-rdata(1,startdate:enddate)).^2,2);
    stdsimu= (1/T*rresult(end))^0.5;
 %% Simulations
% 
% 
% rt=zeros(nsample,T); %store simulation result
% rsimu=zeros(nsimu,T); %store path simulation results
% rsimu(:,1)= rdata(startdate) ; %input initial value
% 
% gtbar=zeros(nsimu,T); %store posterior means over time
% gtbar(:,1)=rdata(startdate);
% sigmatbar=zeros(nsimu,T); %store posterior volatility over time
% sigmag=sigmatbar(:,1);
% stdsimu=zeros(1,nsample);
% 
% for k=1:nsample
%     sigmag=sigmagsample(1,k);
%     n1 = randn(nsimu,T);
%     n2 = randn(nsimu, T);
%     for t = 2:T
%         A=cumsum(rsimu,2);
%         gtbar(:,t)=1/(t+sigma^2/sigmag^2-1)*A(:,t);
%         sigmatbar(:,t)=1/(1/sigmag^2+t/sigma^2);
%         rsimu(:,t)=sigma*n1(:,t)+gtbar(:,t)+sigmatbar(:,t).*n2(:,t);
%         %rsimu(:,t)=rsimu(:,t-1)+sigma*n1(:,t)+sigmatbar(:,t).*n2(:,t);
%     end
%    B=cumsum(rsimu);
%    rt(k,:)=1/nsimu*B(end,:); 
%    rresult=cumsum((rt(k,:)-rdata(1,startdate:enddate)).^2,2);
%    stdsimu(1,k)= (1/T*rresult(end))^0.5;
% end
%%
% %plot sigmagsample vs std
% figure
% plot(sigmagsample,stdsimu);
% title('standard deviation against sigmag');
% xlabel('sigmag');
% ylabel('standard deviation');
% 
% %find minimum std
% [x,y]=find(stdsimu==min(stdsimu));
% stdmin=sigmagsample(x,y)

%%
% plot rdata vs rsimu
figure
t=(tUS(startdate:enddate))';
plot(t, rdata(startdate:enddate));
datetick('x','yyyy-mm-dd');
title('rdata vs rsimu');
xlabel('Date');
ylabel('r');
hold on;
plot(t, rsimu(:,:), 'r');
hold off;
legend('rdata', 'rsimu');

figure
t=(tUS(startdate:enddate))';
plot(t, sigmatbar);
datetick('x','yyyy-mm-dd');
title('sigmatbar against time');
xlabel('Date');
ylabel('parameters');
% hold on;
% plot(t(2:end), gtbar(2:end), 'r');
% hold off;
legend('sigmatbar');

figure
t=(tUS(startdate:enddate))';
%plot(t, sigmatbar);
plot(t, theta, 'r');
datetick('x','yyyy-mm-dd');
title('theta vs date');
xlabel('Date');
ylabel('parameters');
% hold on;
% hold off;
legend('theta');