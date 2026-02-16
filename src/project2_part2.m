%% Question 2
clc; clear; close all;

s = tf('s');
G = 2 / (s^2 + 2*s + 4); 
H = 2 / (s + 1); % from h(t)=2e^-t graph
k_val = -2;
Kp = 0.5;
C = Kp;

%% a) Closed Loop Transfer Function
% T(s) = C(s) * G(s) * (1 - H(S) * k)/(1 - H(s)*k + H(s) * C(s) * G(s))
% T(s) calculated by hand
num_T = C * G * (1 - H * k_val);
den_T = 1 - H * k_val + H * C * G;

% used minreal to erase extra terms like s^5
T = minreal(num_T / den_T, 1e-4);

fprintf('Closed Loop Transfer Function with k = -2 and Kp = 0.5\n');
tf(T)

%% b) Overshoot and Settling Time
info = stepinfo(T);
[Wn, Zeta, P] = damp(T);

% We find only complex poles that their imaginary part isn't 0
complex_indices = find(imag(P) ~= 0);
if ~isempty(complex_indices)
    % Choose the dominant pole with the lowest frequency
    [~, min_idx] = min(Wn(complex_indices));
    dom_idx = complex_indices(min_idx);
    z_dom = Zeta(dom_idx);
    wn_dom = Wn(dom_idx);
else
    % If there is no complex pole, then choose the first pole
    z_dom = Zeta(1);
    wn_dom = Wn(1);
end

fprintf('\nExact Results (stepinfo):\n');
fprintf('Overshoot: %% %.2f\n', info.Overshoot);
fprintf('Settling Time: %.2f s\n', info.SettlingTime);

fprintf('\nApproximate Results:(Formulas)\n');
OS_approx = 100 * exp((-z_dom * pi) / sqrt(1 - z_dom^2));
Ts_approx = 4 / (z_dom * wn_dom);
fprintf('Dominant Pole Zeta: %.4f\n', z_dom);
fprintf('Approximate Overshoot: %% %.2f\n', OS_approx);
fprintf('Approximate Settling Time: %.2f s\n', Ts_approx);
% Since there is a zero in the system the formula may be deviate

%% c) Steady-State

% Final Value Theorem: s*T(s)*1/s = T(0))
yss = dcgain(T); 
fprintf('\nlim(t->inf) y(t) for k = -2: %.4f\n', yss);

% Find k values to get zero ss error
% T(0) should be 1
syms k_sym
G0 = 0.5; % G(0)
H0 = 2;   % H(0)
C0 = 0.5; % Kp
% (C0*G0*(1 - H0*k_sym)) / (1 - H0*k_sym + H0*C0*G0) = 1
T0_eq = (C0*G0*(1 - H0*k_sym)) / (1 - H0*k_sym + H0*C0*G0);
k_zero_error = double(solve(T0_eq == 1, k_sym));

fprintf('K values to get zero SS error %.4f\n', k_zero_error);

% Stability Test
T_test = minreal((C*G*(1-H*k_zero_error))/(1-H*k_zero_error + H*C*G), 1e-4);
if all(real(pole(T_test)) < 0)
    fprintf('The system is stable for the k value\n');
else
    fprintf('The system is unstable although the system mathematically has zero ss error\n');
end

%% d) Kp Range
k_d = 1;
kp_vec = 1.0:0.001:1.5;
valid_kp = [];

for kp_test = kp_vec
    sys_d = minreal((kp_test*G*(1-H*k_d)) / (1-H*k_d + H*kp_test*G), 1e-4);
    poles_d = pole(sys_d);
    % All poles real part smaller than -0.2?
    if all(real(poles_d) < -0.2)
        valid_kp = [valid_kp, kp_test];
    end
end

if ~isempty(valid_kp)
    fprintf('\nValid Kp Range [%.3f, %.3f]\n', min(valid_kp), max(valid_kp));
    figure; pzmap(minreal((min(valid_kp)*G*(1-H*k_d))/(1-H*k_d+H*min(valid_kp)*G)));
    grid on; title(['Kp = ', num2str(min(valid_kp)), ' Pole-Zero Map']);
else
    disp('No Kp was found');
end

figure; step(T); title('k=-2, Kp=0.5 Step Response'); grid on;