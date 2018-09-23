%% LAB 03 - Linearity

%Marc Dominic Cabote
%CECS 463 Fall 2018

close all; clear all; format compact; clc; %clear init
disp('Lab #3 - Linearity');
str = datestr(now); fprintf('MATLAB time stamp: %s\n', str);
disp(' ');

%% Part 1
disp('Part 1');
%a
disp('(a) 2x(n)+1, x1(n)=n, x2(n) = u(n-2), n=[-10,10], a=b=1');
n = -10:10; a=1; b=1;
x1 = n; x2 = stepseq(2,-10,10);
y = a*x1 + b*x2; 
LHS = 2*y + 1; %original equation
RHS = a*(2*x1+1) + b*(2*x2+1);
diff = LHS - RHS;

if (all(diff)==0)
    disp('System is Linear');
else
    disp('System is Non-Linear');
end
disp(' ');

%b
clear all; %clear for next problem
disp('(b) n(x(n-1)) ,x1(n)=cos(0.1pi*n), x2(n)=(0.8)^n, n=[0,20], a=0.5, b=-0.5');
n = 0:20; a=0.5; b=-0.5;
x1 = cos(0.1*pi*n); x2 = (0.8).^n; 
y = a*x1 + b*x2;
LHS = n.*(y-1);
RHS = a*(2*x1+1) + b*(2*x2+1);
diff = LHS - RHS;

if (all(diff)==0)
    disp('System is Linear');
else
    disp('System is Non-Linear');
end
disp(' ');

%c
clear all;
disp('(c) x(n)*x(n-1), x1(n)=1/n, x2(n)=n, n=[-10,10], a=1, b=2');
n=-10:10; a=1; b=2;
x1=1./n; x2=n;
y=a*x1 + b*x2;
LHS = y.*(y-1);
RHS = a*(2*x1+1) + b*(2*x2+1);
diff = LHS - RHS;

if (all(diff)==0)
    disp('System is Linear');
else
    disp('System is Non-Linear');
end
disp(' ');

%d
clear all;
disp('(d) cos(0.5Pi*x(n)), x1(n)=d(n), x2(n)=u(n+2), n=[-10,10], a=2, b=1');
n=-10:10; a=2; b=1;
x1 = impseq(0,-10,10);
x2 = stepseq(-2,-10,10);
y=a*x1 + b*x2;
LHS = cos(0.5*pi*y);
RHS = a*(2*x1+1) + b*(2*x2+1);
diff = LHS - RHS;

if (all(diff)==0)
    disp('System is Linear');
else
    disp('System is Non-Linear');
end
disp(' ');

%% Part 2
disp('Part 2');
%a
disp('(a) 2x(n)+1, x(n)=n, n=[-10,10], k=1');
n=-10:10; k=1; 
x=n;
LHS = (2.*x.*n + 1)-k;
RHS = sigshift(x,n,k);
diff = LHS - RHS;

if (all(diff)==0)
    disp('System is Shift Invariant');
else
    disp('System is Non-Shift Invariant');
end
disp(' ');

%b
clear all;
disp('(b) nx(n-1), x(n)=cos(0.1Pi*n), n=[0:20], k=2');
n=0:20; k=2;
x=cos(0.1*pi*n);
LHS = (n.*(x-1)-k);
RHS = sigshift(x,n,k);
diff = LHS - RHS;

if (all(diff)==0)
    disp('System is Shift Invariant');
else
    disp('System is Non-Shift Invariant');
end
disp(' ');


%c
clear all;
disp('(c) x(n)x(n-1), x(n)=1/n, n=[-10:10], k=-1');
n=-10:10; k=-1;
x=1./n;
LHS = (x.*(x-1))-k;
RHS = sigshift(x,n,k);
diff = LHS - RHS;

if (all(diff)==0)
    disp('System is Shift Invariant');
else
    disp('System is Non-Shift Invariant');
end
disp(' ');

%d
clear all;
disp('(d) cos(0.5*Pi*x(n)), x(n)=d(n), n=[-5:5], k=-2');
n=-5:5; k=-2;
x=impseq(0,-5,5);
LHS = cos(0.5*pi*x)-k;
RHS = sigshift(x,n,k);
diff = LHS - RHS;

if (all(diff)==0)
    disp('System is Shift Invariant');
else
    disp('System is Non-Shift Invariant');
end
disp(' ');

%% Part 3
clear all;
figure(1); clf(1)

n = 0:20;
x_n= 1-(abs(n-10)/5);
T_x= x_n - sigshift(x_n,n,1);
input = exp(-0.25*n);
y = conv(x_n, input);

subplot (3,1,1);grid on;
stem(x_n); title('Excitation');

subplot (3,1,2);grid on;
stem(input); title('Input Sequence');

subplot (3,1,3);grid on;
stem (y); title('Response');

%% Part 4
clear all;
n = -20:20;
figure(2); clf(2);
w1 = pi/4;
w2 = pi/20;
x1 = exp(1j*w1*n);
x2 = exp(1j*w2*n);
plot(x2); grid on;title('Complex Plane'); xlabel('REAL'); ylabel('IMAG');

figure(3); clf(3); 
% pass through a linear system
y = 2*x1 +1;
stem(y); grid on; title('Linear Plot'); xlabel('x'); ylabel('f(x)');

%% Part 5
clear all;
disp('Part 5');
disp('x(t) = cos(Ot) where O=500Pi, Ts = 2x10^-4;');
O = 500*pi;
T = 2*(10^-4);
w =O*T;  
f = w/(2*pi);

fprintf('angular frequency = %2.4f rad/sample\nlinear frequency = %2.4f cyc/samp\n',w,f);


%% Functions
% sigshift funtion
function [y,n] = sigshift(x,m,n0)
% implements y(n) = x(n-n0)
% -------------------------
% [y,n] = sigshift(x,m,n0)
%
n = m+n0; y = x;
end

%impseq function
function [x,n] = impseq(n0,n1,n2)
% Generates x(n) = delta(n-n0); n1 <= n,n0 <= n2
% ----------------------------------------------
% [x,n] = impseq(n0,n1,n2)
%
if ((n0 < n1) | (n0 > n2) | (n1 > n2))
	error('arguments must satisfy n1 <= n0 <= n2')
end
n = [n1:n2];
%x = [zeros(1,(n0-n1)), 1, zeros(1,(n2-n0))];
x = [(n-n0) == 0];
end

%stepseq function
function [x,n] = stepseq(n0,n1,n2)
% Generates x(n) = u(n-n0); n1 <= n,n0 <= n2
% -----------------------------------------
% [x,n] = stepseq(n0,n1,n2)
%
if ((n0 < n1) || (n0 > n2) || (n1 > n2))
	error('arguments must satisfy n1 <= n0 <= n2');
end
n = n1:n2;
%x = [zeros(1,(n0-n1)), ones(1,(n2-n0+1))];
x = (n-n0) >= 0;
end

