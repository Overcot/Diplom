clc; clear;
alpha = [0.319, 0.409, 0.788, 0.818, 0.818, 0.818, 0.818, 0.818,0.818,0.818, 0.818];
x = 1:11;
figure
hold on;
plot(x, alpha, '-o')
title('alpha(s)')
hold off;

gamma = [0.002728, 0.098784, 0.7585, 2.643219, 5.34322, 8.227, 9.891, 11.075, 12.397, 13.494, 15.016];

figure
hold on;
plot(x, gamma, '-o')
title('gamma(s)')
hold off;
