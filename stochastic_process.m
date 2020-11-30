clear all;
close all;

N0 = 20; %initial size of the population
%b = 0.4; % birth rate
m = 0.1; % death rate
dt = 0.001; %time interval

%% run multiple replicate simulations

tmax = 20; %total time
N_t = tmax/dt; %number of time steps

N_repls = 10000; % number of simulations


Ndata = zeros(N_repls,N_t+1);
Ndata(:,1) = N0; %initialize population for each simulation

for i_repl = 1:N_repls
    N = N0; %current population size (reset for each simulation)
    
    for i_t = 1:N_t % how many time steps
        %Ncellborn = sum(rand(N,1)<b*dt); %count number of cell births
        
        %Ncelldead = sum(rand((N - Ncellborn),1)<m*dt);
        
        Ncelldead = sum(rand(N,1)<m*dt);
        N = N - Ncelldead; %update population size
        
        Ndata(i_repl,i_t+1) = N;
    end       
    
end

t = dt*[0:N_t]; %time vector (t = 0, dt, 2*dt, ...)

% Note: these plots show the behavior of the stochastic process. You do not
% need all of these plots for the homework, but you may find them
% interesting.

% make a plot of a few realizations of the population as a function of time
figure(1); 
%plot(t,Ndata(1:min(10,N_repls),:));
plot(t,Ndata(1,:));
xlabel('t'); ylabel('N(t)');
title('Different realizations give different functions N(t)');

% make a plot of the average population size as a function of time
nbar = mean(Ndata,1); % calculate the average population size based off the data

figure(2);
plot(t,nbar,'linewidth',2); %hold on; % computed population size
%plot(t,N0*exp(b*t),'linewidth',2); % average, using the differential equation
xlabel('t'); ylabel('average population size');
hold on
t1 = 0:0.1:tmax;
avgFunc = N0 * exp(-m*t1);
plot(t1,avgFunc,"o");
hold off

popsize = min(min(Ndata)):N0; %vector of population sizes at end of simulation


%nFunc = N_repls*exp(-b*tmax)*(1-exp(-b*tmax)).^(popsize-1); %P_n(t)
nFunc = zeros(length(popsize),1);
for i = 1:length(popsize)
    nFunc(i,1) = N_repls*nchoosek(N0,N0-popsize(i))*exp(-m*N0*tmax)*(exp(m*tmax)-1).^(N0-popsize(i));
end
%plot comparing P_n(t) with distribution from simulations
figure(3);
hist(Ndata(:,end),popsize);
hold on
plot(popsize,nFunc,'r','linewidth',2)
xlabel('pop size'); ylabel('frequency');


% find the growth rate
% lognbar = log(nbar);
% coeff = polyfit(t,lognbar,1);
% slope = coeff(1);
% disp("slop:");disp(slope);