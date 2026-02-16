%% KON 307E - Project 2
%% QUESTION 3: First-Order System Response and Controller/Filter Tasks
% Output files:
%   Q3_step_vs_euler.png
%   Q3_nyquist.png
%   Q3_noise_filter.png

clear; clc; close all;

%% ------------------------------------------------------------
%% (Given/identified) First-order system parameters
s   = tf('s');
K   = 1.273;
tau = 0.318;

% Transfer function
G = K/(tau*s + 1);

fprintf('Q3: G(s) = %.4f / (%.4f s + 1)\n', K, tau);

%% ------------------------------------------------------------
%% 3.1 (Reference) Analytical time response for step input
% For unit step input u(t)=1: y(t)=K(1-exp(-t/tau))
% This part is mainly for completeness & checking.

t_ana = 0:0.001:5;
y_ana = K*(1 - exp(-t_ana/tau));

y_final = K;
t_63 = tau;  % time constant property

fprintf('Steady-state value (unit step) y(∞) = K = %.4f\n', y_final);
fprintf('Time constant tau = %.4f s (y(tau)=0.632*K)\n', t_63);

%% ------------------------------------------------------------
%% 3.2 Step response comparison: MATLAB step() vs Explicit Euler
% MATLAB built-in step response
[y_step, t_step] = step(G, 5);

% Explicit Euler on state form: dy/dt = (-y + K)/tau  for unit step input
h = 0.01;
t = 0:h:5;
y = zeros(size(t));
for k = 1:length(t)-1
    y(k+1) = y(k) + h*((-y(k) + K)/tau);
end

fig32 = figure('Color','w');
plot(t_step, y_step, 'LineWidth', 1.8); hold on;
plot(t, y, '--', 'LineWidth', 1.8);
grid on;
title('Step response: analytical (step) vs Explicit Euler');
xlabel('Time (s)'); ylabel('Output y(t)');
legend('MATLAB step()', 'Explicit Euler', 'Location', 'southeast');

ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1;
ax.XColor = 'k'; ax.YColor = 'k';

exportgraphics(fig32, 'Q3_step_vs_euler.png', 'Resolution', 300);

%% ------------------------------------------------------------
%% 3.3 Nyquist plot for open-loop L(s) = k G(s), k=1
k_gain = 1;
L = k_gain*G;

fig33 = figure('Color','w');
nyquist(L);
grid on;
title('Nyquist plot of L(s)=kG(s), k=1');
xlabel('Real Axis'); ylabel('Imaginary Axis');



exportgraphics(fig33, 'Q3_nyquist.png', 'Resolution', 300);

%% ------------------------------------------------------------
%% 3.4 Filter design to attenuate 100 Hz noise by factor of 10
% Noise: n(t)=A*sin(2*pi*f*t), f=100 Hz
% First-order LPF: F(s)=1/(T s + 1)
% Magnitude: |F(jw)| = 1/sqrt(1+(wT)^2)
% Requirement: |F(jw_n)| = 0.1  at w_n=2*pi*100
% => 1/sqrt(1+(wT)^2)=0.1  -> 1+(wT)^2=100 -> wT=sqrt(99) -> T=sqrt(99)/w

A = 0.25;
f = 100;
w = 2*pi*f;

T = sqrt(99)/w;         % theoretical value
F = tf(1, [T 1]);        % low-pass filter

fprintf('Designed LPF: F(s)=1/(T s + 1),  T = %.6f s\n', T);

% verify magnitude at 100 Hz
mag_100 = abs(freqresp(F, w));  % |F(jw)|
fprintf('|F(jw)| at 100 Hz = %.4f (target 0.1)\n', mag_100);

% simulate filtering in time domain (use long enough time for steady-state)
dt = 1e-4;
t_end = 2;
t_noise = 0:dt:t_end;

noise = A*sin(w*t_noise);
filtered = lsim(F, noise, t_noise);

% show a steady-state window for clean view
t0 = 1.95;  t1 = 2.00;
idx = (t_noise >= t0) & (t_noise <= t1);

figure('Color','w');
plot(t_noise(idx), noise(idx), 'LineWidth', 1.6); hold on;
plot(t_noise(idx), filtered(idx), '--', 'LineWidth', 1.8);
grid on;
title('Noise attenuation at 100 Hz (steady-state window)');
xlabel('Time (s)'); ylabel('Amplitude');
legend('Before filter', 'After filter', 'Location', 'northeast');

ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1;
ax.XColor = 'k'; ax.YColor = 'k';

exportgraphics(gcf, 'Q3_noise_filter.png', 'Resolution', 300);

% Optional: steady-state attenuation estimation from last window (rough)
A_before = (max(noise(idx)) - min(noise(idx)))/2;
A_after  = (max(filtered(idx)) - min(filtered(idx)))/2;
att_ss = A_after / A_before;

fprintf('Steady-state amplitude before = %.4f\n', A_before);
fprintf('Steady-state amplitude after  = %.4f\n', A_after);
fprintf('Estimated attenuation ratio (after/before) = %.4f\n', att_ss);

%% ------------------------------------------------------------
disp('Done. Exported: Q3_step_vs_euler.png, Q3_nyquist.png, Q3_noise_filter.png');
