% Info from graph
t_start = 1;         % Input time
t_at_peak = 3;       % Peak time
y_ss = -0.8;         % Steady State
y_peak = -1.0;       % Max peak value
u_amp = 1;           % Input step amplitude

% Calculations of peak time, Overshoot and Gain
tp = t_at_peak - t_start;                  % Peak time
OS_percent = (abs(y_peak - y_ss) / abs(y_ss)) * 100; % Overshoot (%)
K = y_ss / u_amp;                          % Gain

% Calculations of Damping ratio and Natural Frequency
OS_decimal = OS_percent / 100;
zeta = -log(OS_decimal) / sqrt(pi^2 + (log(OS_decimal))^2);
wn = pi / (tp * sqrt(1 - zeta^2));

% Results
fprintf('Overshoot: %%.%.2f\n', OS_percent);
fprintf('Peak Time (tp): %.2f s\n', tp);
fprintf('Gain (K): %.2f\n', K);
fprintf('Damping Ratio (zeta): %.4f\n', zeta);
fprintf('Natural Frequency (wn): %.4f rad/s\n', wn);


% Physical Parameters
k = -1 / K;          % Spring constant
M = k / (wn^2);          % Mass
b = 2 * zeta * wn * M;   % viscous friction coefficient

fprintf('Parameters:\n');
fprintf('M = %.3f kg, b = %.3f Ns/m, k = %.3f N/m\n', M, b, k);

% Open Loop G(s) = -1 / (M*s^2 + b*s + k)
num = -1;
den = [M b k];
G = tf(num, den);

% Impulse Response 
figure;
impulse(G);
grid on;
title('Impulse Response');
xlabel('Time (sn)');
ylabel('Amplitude');

% New parameters to halve peak time
M_new = M / 4;   % Mass
b_new = b / 2;   % viscous friction coefficient
k_new = k;       % spring constant

% Transfer Functions G(s) = -1 / (Ms^2 + bs + k)
num = -1;
G_old = tf(num, [M b k]);
G_new = tf(num, [M_new b_new k_new]);

% Delaying Step Input
t = 0:0.01:10;
u = (t >= 1); % Input Step at t=1

% Calculating new responses
y_old = lsim(G_old, u, t);
y_new = lsim(G_new, u, t);

% Plotting
figure;
plot(t, y_old, 'b', 'LineWidth', 1.5); hold on;
plot(t, y_new, 'r--', 'LineWidth', 1.5);
grid on;
legend('Old System(tp = 2s)', 'New System (tp = 1s)');
title('Halving the peak time');
xlabel('Time (sn)');
ylabel('Amplitude (y)');