%% LAB #0 - The complex plane

%Marc Dominic Cabote
%CECS 463 Fall 2018

close all; clear all; format compact; clc; %clear init
disp('Lab #0 - The Complex Plane');
str = datestr(now); fprintf('MATLAB time stamo: %s\n', str);
disp(' ');

%% Problem 1
disp('Problem 1');
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
disp('(a)');
x=1-2j; xreal=real(x); ximag=imag(x);
y=-1+1j; yreal=real(y); yimag=imag(y);
a=x+y; b=x-y; c=x*y; d=x/y; p=x^y; %quantities to verify

%check if a is correct
aCheck=(xreal+yreal)+1j*(ximag+yimag); diff=a-aCheck;
if(diff==0)
    disp('quantity a is verified');
else
    disp('error in calculation');
end
%check if b is correct
bCheck=(xreal-yreal)+1j*(ximag-yimag); diff=b-bCheck;
if(diff==0)
    disp('quantity b is verified');
else
    disp('error in calculation');
end
%check if c is verified
cCheck=((xreal*yreal)-(ximag*yimag))+1j*((xreal*yimag)+(ximag*yreal)); 
diff=c-cCheck;
if(diff==0)
    disp('quantity c is verified');
else
    disp('error in calculation');
end
%check if d is verified
dCheck=((xreal*yreal)+(ximag*yimag))+1j*((ximag*yreal)-(xreal*yimag));
dCheck=dCheck/((yreal^2)+(yimag^2));
diff=d-dCheck;
if(diff==0)
    disp('quantity d is verified');
else
    disp('error in calculation');
end
%check if 
xmag=sqrt(xreal^2+ximag^2); xang=angle(x);%angle and magnitude of x
pCheck=exp(y*(log(xmag)+1j*(xang))); diff=p-pCheck;
if(diff==0)
    disp('quantity p is verified');
else
    disp('error in calculation');
end
disp(' ');
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
figure(1); clf(1);
hold on; grid on;
title('Individual points from (a) with Vectors');%graoh title
xlabel('REAL'); ylabel('IMAGINARY');%axislabels
axis([-4,4,-4,4]); %set center
plot((0),(0),'k+');%origin
plot(real(x),imag(x),'Marker','*');
plot(real(y),imag(y),'Marker','*');
plot(real(a),imag(a),'Marker','*');
plot(real(b),imag(b),'Marker','*');
plot(real(c),imag(c),'Marker','*');
plot(real(d),imag(d),'Marker','*');
plot(real(p),imag(p),'Marker','*');
legend ('origin','x','y', 'a', 'b', 'c', 'b', 'p');
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
disp('(c)');
mag=abs(a); ang=angle(a)*(180/pi);
fprintf('Magnitude of a is %4.2f with angle of %4.2f degrees\n',mag,ang);
mag=abs(b); ang=angle(b)*(180/pi);
fprintf('Magnitude of b is %4.2f with angle of %4.2f degrees\n',mag,ang);
mag=abs(c); ang=angle(c)*(180/pi);
fprintf('Magnitude of c is %4.2f with angle of %4.2f degrees\n',mag,ang);
mag=abs(d); ang=angle(d)*(180/pi);
fprintf('Magnitude of a is %4.2f with angle of %4.2f degrees\n',mag,ang);
mag=abs(p); ang=angle(p)*(180/pi);
fprintf('Magnitude of a is %4.2f with angle of %4.2f degrees\n',mag,ang);
disp(' ');

plot([0,real(x)],[0,imag(x)]);
plot([0,real(y)],[0,imag(y)]);
plot([0,real(a)],[0,imag(a)]);
plot([0,real(b)],[0,imag(b)]);
plot([0,real(c)],[0,imag(c)]);
plot([0,real(d)],[0,imag(d)]);
plot([0,real(p)],[0,imag(p)]);
legend ('origin','x','y', 'a', 'b', 'c', 'b', 'p');
hold off;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%create four random points each quadrant
disp('(d)');
figure(2);clf(2);
hold on; grid on;

q1=+randi(5)+1j*(randi(5)); %Q1
q2=-randi(5)+1j*(randi(5)); %Q2
q3=-randi(5)-1j*(randi(5)); %Q3
q4=+randi(5)-1j*(randi(5)); %Q4

fprintf('Q1=%4.2f+%4.2fj\n', real(q1),imag(q1));
fprintf('Q2=%4.2f+%4.2fj\n', real(q2),imag(q2));
fprintf('Q3=%4.2f+%4.2fj\n', real(q3),imag(q3));
fprintf('Q4=%4.2f+%4.2fj\n', real(q4),imag(q4));

area_triangle=@(x1,x2,x3) 0.5*abs(real(x1)*(imag(x2)-imag(x3)) ...
+ real(x2)*(imag(x3)-imag(x1)) + real(x3)*(imag(x1)-imag(x2)));
area_quad=@(x1,x2,x3,x4) area_triangle(x1,x2,x3) + area_triangle(x1,x3,x4);
Area = area_quad(q1,q2,q3,q4);

title(sprintf('AREA =%4.2f',Area));%graph title which is area
xlabel('REAL'); ylabel('IMAGINARY');%axislabels
axis([-5,5,-5,5]); %set center
plot((0),(0),'k+');%origin

plot(real(q1),imag(q1),'Marker','*');
plot(real(q2),imag(q2),'Marker','*');
plot(real(q3),imag(q3),'Marker','*');
plot(real(q4),imag(q4),'Marker','*');

plot([real(q1),real(q2)],[imag(q1),imag(q2)]);
plot([real(q2),real(q3)],[imag(q2),imag(q3)]);
plot([real(q3),real(q4)],[imag(q3),imag(q4)]);
plot([real(q4),real(q1)],[imag(q4),imag(q1)]);

%% Problem 2
figure(3); clf(3); grid on; 

title('Sinusoids as phasors');%graoh title
xlabel('REAL'); ylabel('IMAGINARY');%axislabels
axis ([-4,4,-4,4]);

x=3*exp(1j*45*(pi/180)); xreal=real(x); ximag=imag(x);
y=2*exp(1j*(-150-90)*(pi/180)); yreal=real(y); yimag=imag(y);

hold on;
plot(0,0,'+');%origin
plot(real(x),imag(x),'b*');
plot(real(y),imag(y),'r*');
plot(real(x+y),imag(x+y),'g*')
plot([0,xreal],[0,ximag],'b');
plot([0,yreal],[0,yimag],'r');
disp(' ');
%% Problem 3
disp('Problem 3');
fprintf('Magnitude(x+y)= %4.2f phase_angle = %4.2f',abs(x+y),angle(x+y)*(180/pi));
disp(' degrees');
plot([0,xreal+yreal],[0,ximag+yimag],'g');
legend ('origin','x(t)','y(t)','x(t)+y(t)');
disp(' ');
hold off;

%% Problem 4
disp('Problem 4'); 
Offset=(-2.75+0.75)/2;
A=(0.75-(-2.75))/2;%amplitude
angFreq=(2*pi)/(0.04);%angular frequency
phase=(0-0.018)*(360/0.04);

fprintf('Offset=%4.2f;\nAmplitude=%4.2f;\n',Offset,A);
fprintf('Angular Frequency=%4.2f;\nphase=%4.2f;\n',angFreq,phase);




