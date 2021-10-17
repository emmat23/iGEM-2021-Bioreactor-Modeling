%Modeling the bacterial growth in the bioreactor

%Growth Kinetics:

%N = Ni*(e^(u*t))
%u = specific growth rate of organism
%t = time
%N = amount of bacteria over time t (number of cells over time)
%Ni = initial amount of bacteria (number of cells)

%Growth rate found using doubling time
%td = ln(2)/u
%deniz said doubling time is prob 20 min

%Modeling reactions in bioreactor:

%Rate equation for cell mass production in reactor
%rate = (u*x) - ((f/V)*x)
%u = specific growth rate of organism
%x = cell mass (g/L)
%f = volume flow (L/s)
%V = volume (L)

td = 20; %doubling time is 20 minutes
u = findGrowth(td);

Ni = 100; %number of cells
t_n = linspace(0,60,85); %minutes
N = ecoliGrowth(u,t_n,Ni); %number of cells at time t


V = 1.5; %volume (liters)
fmin = (170/1000); %flow rate of motor(L/min)
fmax = (460/1000); %volume flow in L/min
x0 = 100; %number of initial cells
tspan = [0, 60]; %minutes

[t,x] = ode45(@(t,x) monodModel(x,fmax,u),tspan,x0);

dNdt = monodModel(x,fmax,u);

n = length(dNdt);

for i = 1:n
    K = carry_cap(dNdt(i),u,N(i));
end
% 
%Plot Growth Rate

figure
plot(t, N,'b-')
ylim([0 900]);
hold on
yline(K,'r-', 'Carrying Capacity: 800')
hold off
title('Growth Rate of E.Coli')
xlabel('Time (min)')
ylabel('Number of Cells of Bacteria')

%Subfunction for the MonodModel
function dNdt = monodModel(x,f,u)
    V = 1.5;
    dNdt = (u*x) - ((f/V)*x);
end

%Subfunction to find the specific growth rate
function GR = findGrowth(td)
    GR = (log(2))/td;
end

%Subfunction to find the growth of E.Coli
function ecoli = ecoliGrowth(u,t,Ni)
    ecoli = Ni*(exp(u*t));
end

%Subfunction for carrying capacity
function K = carry_cap(dNdt,u,N)
    K = N/(1-((dNdt)/(u*N)));
end