clc; clear;
alpha = [1, 0.95, 0.85, 0.8, 0.7, 0.5, 0.3, 0];
x = 1:8;
figure
hold on;
plot(x, alpha, '-o')
title('alpha(s) from diplom')
hold off;

gamma(x) = ((x-1).*(7-(x-1)))/(8^2)
figure
hold on;
plot(x, gamma, '-o')
title('gamma(s) from diplom')
hold off;

x0 = [280026, 118907, 33972, 10172, 2456, 993, 483, 3];
figure
hold on;
plot(x, x0, '-o')
title('x0(s) from diplom')
hold off;

alpha = [0.319, 0.409, 0.788, 0.818, 0.818, 0.818, 0.818, 0.818,0.818,0.818, 0.818];
x=1:11;
figure
hold on;
plot(x, alpha, '-o')
title('alpha(s) from new data')
hold off;

gamma = [0.002728, 0.098784, 0.7585, 2.643219, 5.34322, 8.227, 9.891, 11.075, 12.397, 13.494, 15.016];
figure
hold on;
plot(x, gamma, '-o')
title('gamma(s) from new data')
hold off;