# MAE384_FinalProject
% Parameters for simulation
h = 1; % Time step
T = 100; % Total time
t = 0:h:T; % Time vector
N = 1000; % Total population
initial_conditions = [990, 10, 0]; % [S(0), I(0), R(0)]
parameters = [0.3, 0.1; 1, 0.1;2, 0.2];%i. Seasonal Influenza: β = 0.3, γ = 0.1 ii. COVID-19: β = 1, γ = 0.1 iii. Measles: β = 2, γ = 0.2
dSdt= @(S,I) -B/N*S*I;
dIdt= @(S,I) B/N*S*I-y*I;
dRdt= @(I) y*I;


% Newton Linear Interpolation
f(x)=f(x0)+(f(x1)-f(x0))/(x1-x0)*(x-x0);
%Newton Quadratic Interpolation
f(x)=f(x0)+(f(x1)-f(x0))/(x1-x0)*(x-x0)+((f(x2)-f(x1))/(x2-x1)-(f(x1)-f(x0))/(x1-x0))/(x2-x0)*(x-x0)*(x-x1);

%% Part 3:Least Squares
I0_true=10
B_true=.3
% 30 days
x=t(1:30);
y=log(II(1:30));
n=30;
N=1000;
S0=990;
gamma=0.1;
a1=(n*sum(x.*y)-sum(x)*sum(y))/(n*sum(x.^2)-(sum(x))^2);
a0=sum(y)/n-a1*sum(x)/n;
k=a1;
I0_30=exp(a0)
B_30=(k+gamma)*N/S0
% 10 days 
x=t(1:10);
y=log(II(1:10));
n=10;
a1=(n*sum(x.*y)-sum(x)*sum(y))/(n*sum(x.^2)-(sum(x))^2);
a0=sum(y)/n-a1*sum(x)/n;
k=a1;
I0_10=exp(a0)
B_10=(k+gamma)*N/S0
