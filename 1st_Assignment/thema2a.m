clc, clearvars, close all;

%Generating samples for q and u

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
dt = 1e-4;
tspan = 0:dt:20;
[~, x] = ode45(f, tspan, x0);

Ts = 0.1;
t_samples = 0:Ts:20;

u_samples = u(t_samples');
q_samples = interp1(tspan, x(:,1), t_samples');
q_dot_samples = interp1(tspan, x(:,2), t_samples');

%Filtering 

%Lambda = tf([1 1], 1); % s+1 
%H = tf([1 0], [1  1]); % s / (s+1) 
Lambda = tf([1 3 2], 1); % s^2+3s+2
H = tf([1 0], [1 3 2]); % s / (s^2 + 3s + 2)
J1 = lsim(1/Lambda, -q_dot_samples, t_samples');
J2 = lsim(1/Lambda, -q_samples, t_samples');
J3 = lsim(1/Lambda, u_samples, t_samples');
J = [J1 J2 J3];

Y = lsim(H, q_dot_samples, t_samples');

theta = J \ Y;

%From thema1
A = [0 1;
    -theta(2) -theta(1)];
B = [0;
    theta(3)];

L_hat = g / theta(2);
m_hat = 1 / (theta(3) * L_hat^2);
c_hat = theta(1) * m_hat * L_hat^2;

fprintf("Οι εκτιμώμενες τιμές των παραμέτρων είναι:\n");
fprintf("Μάζα m: %.4f\n", m_hat);
fprintf("Μήκος εκκρεμούς L: %.4f\n", L_hat);
fprintf("Συντελεστής απόσβεσης c: %.4f\n", c_hat);

 

%Simulation of model
u = @(t) A0 * sin(w * t);
f = @(t, xhat) A*xhat + B*u(t);
x0 = [0; 0];
[~, xhat] = ode45(f, tspan, x0);
qhat_samples = interp1(tspan, xhat(:, 1), t_samples');

figure;
plot(t_samples, q_samples);  
xlabel('Χρόνος t σε s');
ylabel('Πραγματική γωνία εκτροπής q(t) σε rad');
title('Πραγματική απόκριση γωνίας εκτροπής q(t)');
grid on;

figure;
plot(t_samples, qhat_samples);
xlabel('Χρόνος t σε s');
ylabel('Εκτιμώμενη γωνία εκτροπής qhat(t) σε rad');
title("Εκτιμώμενη απόκριση γωνίας εκτροπής qhat(t)");
grid on;

e = q_samples - qhat_samples;
figure;
plot(t_samples, e);
xlabel('Χρόνος t σε s');
ylabel('Σφάλμα e(t) σε rad');
title("Τιμές σφάλματος e(t) = q(t) - qhat(t)");
grid on;