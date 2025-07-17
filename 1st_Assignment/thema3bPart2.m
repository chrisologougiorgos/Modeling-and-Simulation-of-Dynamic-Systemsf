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

Ts = [0.01:0.01:0.5];
parameters_error = zeros(1, length(Ts));

for i=1:length(Ts)
    
    t_samples = 0:Ts(i):20;
    u_samples = u(t_samples');
    q_samples = interp1(tspan, x(:,1), t_samples');
    
    %Filtering 
    Lambda = tf([1 3 2], 1); % s^2 + 3s + 2
    H1 = tf([1 0], [1 3 2]); % s / (s^2 + 3s + 2)
    H2 = tf([1 0 0], [1 3 2]); % s^2 / (s^2 + 3s + 2)
    J1 = lsim(H1, -q_samples, t_samples');
    J2 = lsim(1/Lambda, -q_samples, t_samples');
    J3 = lsim(1/Lambda, u_samples, t_samples');
    J = [J1 J2 J3];
    
    Y = lsim(H2, q_samples, t_samples');
   
    theta = J \ Y;
    
    %From thema1
    L_hat = g / theta(2);
    m_hat = 1 / (theta(3) * L_hat^2);
    c_hat = theta(1) * m_hat * L_hat^2;
    
    fprintf("Τιμή περιόδου δειγματοληψίας Ts: %f\n", Ts(i))
    fprintf("Οι εκτιμώμενες τιμές των παραμέτρων είναι:\n");
    fprintf("Μάζα m: %.4f\n", m_hat);
    fprintf("Μήκος εκκρεμούς L: %.4f\n", L_hat);
    fprintf("Συντελεστής απόσβεσης c: %.4f\n\n", c_hat);

    parameters_error(i) = (m - m_hat)^2 + (L - L_hat)^2 + (c - c_hat)^2;

end

figure;
plot(Ts, parameters_error);
xlabel("Περίοδος δειγματοληψίας Ts σε s");
ylabel("Tετραγωνικό σφάλμα εκτίμησης παραμέτρων");
title("Σφάλμα εκτίμησης παραμέτρων ανά τιμή Ts");
grid on;