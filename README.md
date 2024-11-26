# MAE384_FinalProject
% Parameters for simulation
h = 1; % Time step
T = 100; % Total time
t = 0:h:T; % Time vector
N = 1000; % Total population
initial_conditions = [990, 10, 0]; % [S(0), I(0), R(0)]
parameters = [0.3, 0.1; 1, 0.1;2, 0.2];%i. Seasonal Influenza: β = 0.3, γ = 0.1 ii. COVID-19: β = 1, γ = 0.1 iii. Measles: β = 2, γ = 0.2
dSdt= @(B,S,I,y) -B/N*S*I;
dIdt= @(B,S,I,y) B/N*S*I-y*I;
dRdt= @(B,S,I,y) y*I;


% Newton Linear Interpolation
f(x)=f(x0)+(f(x1)-f(x0))/(x1-x0)*(x-x0);
%Newton Quadratic Interpolation
f(x)=f(x0)+(f(x1)-f(x0))/(x1-x0)*(x-x0)+((f(x2)-f(x1))/(x2-x1)-(f(x1)-f(x0))/(x1-x0))/(x2-x0)*(x-x0)*(x-x1);
