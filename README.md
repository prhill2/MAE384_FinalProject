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
% 5 days 
x=t(1:5);
y=log(II(1:5));
n=10;
a1=(n*sum(x.*y)-sum(x)*sum(y))/(n*sum(x.^2)-(sum(x))^2);
a0=sum(y)/n-a1*sum(x)/n;
k=a1;
I0_5=exp(a0)
B_5=(k+gamma)*N/S0


%% Part 1 Practice Code

% Parameters
N = 1000; % Total population
h = 1;    % Step size
T = 100;  % Total time
N_t = T / h; % Number of time steps


% Initial conditions
S(1) = 990; 
I(1) = 10;  
R(1) = 0;   

% Gamma and Beta for Influenza
B = 0.3;  
y = 0.1; 


% Runge-Kutta

for i = 1:N_t
    
    dSdt = @(S, I) -B / N * S * I;
    dIdt = @(S, I) B / N * S * I - y * I;
    dRdt = @(I) y * I;

    
    k_1S = dSdt(S(i), I(i));
    k_2S = dSdt(S(i) + 0.5 * h , I(i) + 0.5 * h * k_1S);
    k_3S = dSdt(S(i) + 0.5 * h , I(i) + 0.5 * h * k_2S);
    k_4S = dSdt(S(i) + h , I(i) + h * k_3S);

    
    S(i + 1) = S(i) + (h / 6) * (k_1S + 2 * k_2S + 2 * k_3S + k_4S);


    k_1I = dIdt(S(i), I(i));
    k_2I = dIdt(S(i) + 0.5 * h , I(i) + 0.5 * h * k_1I);
    k_3I = dIdt(S(i) + 0.5 * h , I(i) + 0.5 * h * k_2I);
    k_4I = dIdt(S(i) + h , I(i) + h * k_3I);

    
    I(i + 1) = I(i) + (h / 6) * (k_1I + 2 * k_2I + 2 * k_3I + k_4I);

    
    k_1R = dRdt(I(i));
    k_2R = dRdt(I(i) + 0.5 * h * k_1I);
    k_3R = dRdt(I(i) + 0.5 * h * k_2I);
    k_4R = dRdt(I(i) + h * k_3I);

    
    R(i + 1) = R(i) + (h / 6) * (k_1R + 2 * k_2R + 2 * k_3R + k_4R);

    
    t(i + 1) = t(i) + h; % update t
end

% Plot for Seasonal Influenza
figure;
plot(t, S, 'b*--', t, I, 'r*--', t, R, 'g*--');
title('Seasonal Influenza');
xlabel('Time');
ylabel('Population');
legend('S(t)', 'I(t)', 'R(t)');
grid on;

%%
% Gamma and Beta for Covid-19
B = 1;  
y = 0.1; 


% Runge-Kutta

for i = 1:N_t
    
    dSdt = @(S, I) -B / N * S * I;
    dIdt = @(S, I) B / N * S * I - y * I;
    dRdt = @(I) y * I;

    
    k_1S = dSdt(S(i), I(i));
    k_2S = dSdt(S(i) + 0.5 * h , I(i) + 0.5 * h * k_1S);
    k_3S = dSdt(S(i) + 0.5 * h , I(i) + 0.5 * h * k_2S);
    k_4S = dSdt(S(i) + h , I(i) + h * k_3S);

    
    S(i + 1) = S(i) + (h / 6) * (k_1S + 2 * k_2S + 2 * k_3S + k_4S);


    k_1I = dIdt(S(i), I(i));
    k_2I = dIdt(S(i) + 0.5 * h , I(i) + 0.5 * h * k_1I);
    k_3I = dIdt(S(i) + 0.5 * h , I(i) + 0.5 * h * k_2I);
    k_4I = dIdt(S(i) + h , I(i) + h * k_3I);

    
    I(i + 1) = I(i) + (h / 6) * (k_1I + 2 * k_2I + 2 * k_3I + k_4I);

    
    k_1R = dRdt(I(i));
    k_2R = dRdt(I(i) + 0.5 * h * k_1I);
    k_3R = dRdt(I(i) + 0.5 * h * k_2I);
    k_4R = dRdt(I(i) + h * k_3I);

    
    R(i + 1) = R(i) + (h / 6) * (k_1R + 2 * k_2R + 2 * k_3R + k_4R);

    
    t(i + 1) = t(i) + h; % update t
end

% Plot for Covid-19
figure;
plot(t, S, 'b*--', t, I, 'r*--', t, R, 'g*--');
title('Covid-19');
xlabel('Time');
ylabel('Population');
legend('S(t)', 'I(t)', 'R(t)');
grid on;

%%
% Gamma and Beta for Measles
B = 2;  
y = 0.2; 


% Runge-Kutta

for i = 1:N_t
    
    dSdt = @(S, I) -B / N * S * I;
    dIdt = @(S, I) B / N * S * I - y * I;
    dRdt = @(I) y * I;

    
    k_1S = dSdt(S(i), I(i));
    k_2S = dSdt(S(i) + 0.5 * h , I(i) + 0.5 * h * k_1S);
    k_3S = dSdt(S(i) + 0.5 * h , I(i) + 0.5 * h * k_2S);
    k_4S = dSdt(S(i) + h , I(i) + h * k_3S);

    
    S(i + 1) = S(i) + (h / 6) * (k_1S + 2 * k_2S + 2 * k_3S + k_4S);


    k_1I = dIdt(S(i), I(i));
    k_2I = dIdt(S(i) + 0.5 * h , I(i) + 0.5 * h * k_1I);
    k_3I = dIdt(S(i) + 0.5 * h , I(i) + 0.5 * h * k_2I);
    k_4I = dIdt(S(i) + h , I(i) + h * k_3I);

    
    I(i + 1) = I(i) + (h / 6) * (k_1I + 2 * k_2I + 2 * k_3I + k_4I);

    
    k_1R = dRdt(I(i));
    k_2R = dRdt(I(i) + 0.5 * h * k_1I);
    k_3R = dRdt(I(i) + 0.5 * h * k_2I);
    k_4R = dRdt(I(i) + h * k_3I);

    
    R(i + 1) = R(i) + (h / 6) * (k_1R + 2 * k_2R + 2 * k_3R + k_4R);

    
    t(i + 1) = t(i) + h; % update t
end

% Plot for Measles
figure;
plot(t, S, 'b*--', t, I, 'r*--', t, R, 'g*--');
title('Measles');
xlabel('Time');
ylabel('Population');
legend('S(t)', 'I(t)', 'R(t)');
grid on;
