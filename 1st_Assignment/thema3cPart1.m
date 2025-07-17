clc, clearvars, close all;

%Generating samples for q, q' and u

m = 0.75;
L = 1.25;
c = 0.15;
g = 9.81;
w = 2;

A = [0 1;
    -g/L -c/(m*L^2)];

B = [0 ;
    1/(m*L^2)];

A0 = 1:1:50;
x0 = [0; 0];
dt = 1e-4;
tspan = 0:dt:20;
Ts = 0.1;
t_samples = 0:Ts:20;

parameters_error = zeros(1, length(A0));

for i=1:length(A0)
    
    u = @(t) A0(i) * sin(w * t);
    f = @(t, x) A*x + B*u(t);
    [t, x] = ode45(f, tspan, x0);
    
    u_samples = u(t_samples');
    q_samples = interp1(tspan, x(:,1), t_samples');
    q_dot_samples = interp1(tspan, x(:,2), t_samples');
    
    
    %Filtering 
    %Lambda = tf([1 1], 1); % s+1 
    %H = tf([1 0], [1  1]); % s / (s+1) 
    Lambda = tf([1 3 2], 1); % s^2+3s+2
    H = tf([1 0], [1 3 2]); % s / (s^2 + 3s + 2)
    J1 = lsim(1/Lambda, -q_dot_samples, t_samples);
    J2 = lsim(1/Lambda, -q_samples, t_samples);
    J3 = lsim(1/Lambda, u_samples, t_samples);
    J = [J1 J2 J3];
    
    Y = lsim(H, q_dot_samples, t_samples);
    
    theta = J \ Y;
    
    %From thema1
    L_hat = g / theta(2);
    m_hat = 1 / (theta(3) * L_hat^2);
    c_hat = theta(1) * m_hat * L_hat^2;

    fprintf("Τιμή πλάτους Α: %d\n", A0(i))
    fprintf("Οι εκτιμώμενες τιμές των παραμέτρων είναι:\n");
    fprintf("Μάζα m: %.7f\n", m_hat);
    fprintf("Μήκος εκκρεμούς L: %.7f\n", L_hat);
    fprintf("Συντελεστής απόσβεσης c: %.7f\n\n", c_hat);
    
    parameters_error(i) = (m - m_hat)^2 + (L - L_hat)^2 + (c - c_hat)^2;
end

figure;
plot(A0, parameters_error);
xlabel("Τιμή πλάτους A0 ");
ylabel("Tετραγωνικό σφάλμα εκτίμησης παραμέτρων");
title("Σφάλμα εκτίμησης παραμέτρων ανά τιμή A0");
grid on;