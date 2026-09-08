%% BENDS synthetic example
% This example uses the same made-up segment as the Python and R examples.
% The input signal is used directly. BENDS performs no preprocessing.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
dataPath = fullfile(repoRoot, 'example_data', 'synthetic_segment.csv');
addpath(fullfile(repoRoot, 'MATLAB'));

D = readtable(dataPath);
t = D.t_min;
y = D.signal;

% BENDS settings for this example.
t_N = 35;      % documented event time [min]
t_D = 20;      % maximum timing uncertainty [min]
PreL = 15;     % pre-candidate window length [min]
PostL = 15;    % post-candidate window length [min]
tau = 0.05;    % residual-fluctuation tolerance
q = +1;        % expected shift toward a more positive slope

% Known slope-change time used only to evaluate the synthetic example.
c_true = 30;

[c_star, m1, m2, b1, b2] = BENDS( ...
    t, y, t_N, t_D, PreL, PostL, tau, q);

fprintf('c_star = %.6f min\n', c_star);
fprintf('m1     = %.6f signal units/min\n', m1);
fprintf('m2     = %.6f signal units/min\n', m2);
fprintf('b1     = %.6f signal units\n', b1);
fprintf('b2     = %.6f signal units\n', b2);

% Fitted lines at the selected candidate.
t_pre = t(t >= c_star - PreL & t <= c_star);
x_pre = t_pre - (c_star - PreL);
y_pre_fit = b1 + m1 * x_pre;

t_post = t(t >= c_star & t <= c_star + PostL);
x_post = t_post - c_star;
y_post_fit = b2 + m2 * x_post;

figure('Color', 'w');
plot(t, y, 'k-', 'LineWidth', 1.2);
hold on;
plot(t_pre, y_pre_fit, 'LineWidth', 2);
plot(t_post, y_post_fit, 'LineWidth', 2);
xline(t_N, '--', 'Documented time');
xline(c_true, ':', 'Known synthetic change','LabelVerticalAlignment','bottom');
xline(c_star, '-.', 'BENDS inferred time');
xlabel('Time [min]');
ylabel('Signal');
title('BENDS synthetic example');
grid on;
legend('Signal', 'Pre-candidate fit', 'Post-candidate fit', ...
    'Location', 'best');
