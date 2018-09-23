%% LAB 02 - Discrete Time Signals

%Marc Dominic Cabote
%CECS 463 Fall 2018

close all; clear all; format compact; clc; %clear init
disp('Lab #2 - Discrete Time Signals');
str = datestr(now); fprintf('MATLAB time stamp: %s\n', str);
disp(' ');

%% Part 1
disp('(1)');
n = -3:7;
x=n*0;x(n==0)=2; x(n==2)=1; x(n==3)=-1; x(n==4)=3; 

figure(1);clf(1);
stem(n,x);grid on;

%% Part 2 and 3
disp('(2) & (3)');
n = -3:7;
x=n*0;x(n==0)=2; x(n==2)=1; x(n==3)=-1; x(n==4)=3; 

figure(2);
%y1[n] = x[n-2]
[y1 ,n1] = sigshift (x, n, 2);
subplot(4,1,1);grid on;
stem(n1,y1); title('y1[n] = x[n-2]');

%y2[n] = x[n+1]
[y2, n2] = sigshift (x, n, -1);
subplot(4,1,2);grid on;
stem(n2,y2); title('y2[n] = x[n+1]');

%y3[n] = x[-n]
[y3, n3] = sigshift (x, -n, 0);
subplot(4,1,3);grid on;
stem(n3,y3); title('y3[n] = x[-n]');

%y4[n] = x[-n]
[y4, n4] = sigshift (x, -n, -1);
subplot(4,1,4);grid on;
stem(n4,y4); title('y3[n] = x[-n+1]');


function [y,m] = sigshift(x,n,k)

% implements y(n) = x(n-n0)
% -------------------------
% [y,n] = sigshift(x,m,n0)
m = n+k; 
y = x;
end