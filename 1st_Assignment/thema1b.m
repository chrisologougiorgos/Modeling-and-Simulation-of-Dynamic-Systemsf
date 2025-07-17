clc, clearvars, close all;

m = 0.75;
L = 1.25;
c = 0.15;
g = 9.81;
A0 = 4;
w = 2;

A = [0 1;
    -g/L -c/(m*L^2)];

B = [0 ;
    1/(m*L^2)];

u = @(t) A0 * sin(w * t);
f = @(t, x) A*x + B*u(t);
x0 = [0; 0];
tspan = 0:1e-4:20;
[t, x] = ode45(f, tspan, x0);

figure;
plot(t, x(:,1));  
xlabel('Χρόνος t σε s');
ylabel('Γωνία εκτροπής q(t) σε rad');
title('Απόκριση γωνίας q(t)');
grid on;

figure;
plot(t, x(:,2));
xlabel('Χρόνος t σε s');
ylabel("Γωνιακή ταχύτητα q'(t) σε rad");
title("Απόκριση γωνιακής ταχύτητας q'(t)");
grid on;