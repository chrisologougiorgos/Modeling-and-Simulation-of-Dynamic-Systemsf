clc, clearvars, close all;

%% Generating samples for q, q' and u

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

std = 0.05;
noise = std * randn(length(u_samples), 2);
q_samples_with_noise = q_samples + noise(:, 1);
q_dot_samples_with_noise = q_dot_samples + noise(:, 2);

%% With measurable qdot and first-order filter
%Filtering 
Lambda = tf([1 1], 1); % s+1
H = tf([1 0], [1 1]); % s / (s+1)
J1 = lsim(1/Lambda, -q_dot_samples_with_noise, t_samples');
J2 = lsim(1/Lambda, -q_samples_with_noise, t_samples');
J3 = lsim(1/Lambda, u_samples, t_samples');
J = [J1 J2 J3];

Y = lsim(H, q_dot_samples_with_noise, t_samples');

theta = J \ Y;

%From thema1
A = [0 1;
    -theta(2) -theta(1)];
B = [0;
    theta(3)];

L_hat = g / theta(2);
m_hat = 1 / (theta(3) * L_hat^2);
c_hat = theta(1) * m_hat * L_hat^2;
parameters_error = (m - m_hat)^2 + (L - L_hat)^2 + (c - c_hat)^2;

fprintf("Μοντέλο με μετρήσιμο q_dot και φίλτρο 1ης τάξης");
fprintf("\nΟι εκτιμώμενες τιμές των παραμέτρων είναι:\n");
fprintf("Μάζα m: %.4f\n", m_hat);
fprintf("Μήκος εκκρεμούς L: %.4f\n", L_hat);
fprintf("Συντελεστής απόσβεσης c: %.4f\n", c_hat);
fprintf("Τετραγωνικό σφάλμα εκτίμησης παραμέτρων: %.4f\n", parameters_error)


%% With measurable qdot and second-order filter
%Filtering 

Lambda = tf([1 3 2], 1); % s^2 + 3s + 2
H = tf([1 0], [1 3 2]); % s / (s^2 + 3s + 2)
J1 = lsim(1/Lambda, -q_dot_samples_with_noise, t_samples');
J2 = lsim(1/Lambda, -q_samples_with_noise, t_samples');
J3 = lsim(1/Lambda, u_samples, t_samples');
J = [J1 J2 J3];

Y = lsim(H, q_dot_samples_with_noise, t_samples');

theta = J \ Y;

%From thema1
A = [0 1;
    -theta(2) -theta(1)];
B = [0;
    theta(3)];

L_hat = g / theta(2);
m_hat = 1 / (theta(3) * L_hat^2);
c_hat = theta(1) * m_hat * L_hat^2;
parameters_error = (m - m_hat)^2 + (L - L_hat)^2 + (c - c_hat)^2;

fprintf("\nΜοντέλο με μετρήσιμο q_dot και φίλτρο 2ης τάξης");
fprintf("\nΟι εκτιμώμενες τιμές των παραμέτρων είναι:\n");
fprintf("Μάζα m: %.4f\n", m_hat);
fprintf("Μήκος εκκρεμούς L: %.4f\n", L_hat);
fprintf("Συντελεστής απόσβεσης c: %.4f\n", c_hat);
fprintf("Τετραγωνικό σφάλμα εκτίμησης παραμέτρων: %.4f\n", parameters_error)


%% Without measurable q_dot (same second-order filter as above)

%Filtering 
Lambda = tf([1 3 2], 1); % s^2 + 3s + 2
H1 = tf([1 0], [1 3 2]); % s / (s^2 + 3s + 2)
H2 = tf([1 0 0], [1 3 2]); % s^2 / (s^2 + 3s + 2)
J1 = lsim(H1, -q_samples_with_noise, t_samples');
J2 = lsim(1/Lambda, -q_samples_with_noise, t_samples');
J3 = lsim(1/Lambda, u_samples, t_samples');
J = [J1 J2 J3];

Y = lsim(H2, q_samples_with_noise, t_samples');

theta = J \ Y;

%From thema1
A = [0 1;
    -theta(2) -theta(1)];
B = [0;
    theta(3)];

L_hat = g / theta(2);
m_hat = 1 / (theta(3) * L_hat^2);
c_hat = theta(1) * m_hat * L_hat^2;
parameters_error = (m - m_hat)^2 + (L - L_hat)^2 + (c - c_hat)^2;

fprintf("\nΜοντέλο με μη μετρήσιμο q_dot");
fprintf("\nΟι εκτιμώμενες τιμές των παραμέτρων είναι:\n");
fprintf("Μάζα m: %.4f\n", m_hat);
fprintf("Μήκος εκκρεμούς L: %.4f\n", L_hat);
fprintf("Συντελεστής απόσβεσης c: %.4f\n", c_hat);
fprintf("Τετραγωνικό σφάλμα εκτίμησης παραμέτρων: %.4f\n", parameters_error)
 

