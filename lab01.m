%% LAB 01 - Quadrilaterals

%Marc Dominic Cabote
%CECS 463 Fall 2018

close all; clear all; format compact; clc; %clear init
disp('Lab #1 - Quadrilaterals');
str = datestr(now); fprintf('MATLAB time stamp: %s\n', str);
disp(' ');

%% Quadrilaterals 
disp('Quadrilaterals');
figure(3);clf(3);
hold on; grid on;
axis([-10,10,-10,10]); %set center
plot((0),(0),'k+');%origin

%functions to compute area of quadrilateral
area_triangle=@(x1,x2,x3) 0.5*abs(real(x1)*(imag(x2)-imag(x3)) ...
+ real(x2)*(imag(x3)-imag(x1)) + real(x3)*(imag(x1)-imag(x2)));

area_quad=@(x1,x2,x3,x4) ...
area_triangle(x1,x2,x3) + area_triangle(x1,x3,x4);


%cross product function
xprod = @(x1,x2) real(x1)*imag(x2) - imag(x1)*real(x2);

%generate four random points each quadrant from [0,10]
q1=+randi(10)+1j*(randi(10)); %Q1
q2=-randi(10)+1j*(randi(10)); %Q2
q3=-randi(10)-1j*(randi(10)); %Q3
q4=+randi(10)-1j*(randi(10)); %Q4
%q1=6+1j*4; %test point
%q2=-2+1j*4; %test point
%q3=-2-1j*8; %test point
%q4=2-1j*2; %test point
%generate the vertices from the points
v1=q2-q1;
v2=q3-q2;
v3=q4-q3;
v4=q1-q4;

Area = area_quad(q1,q2,q3,q4);
%find the cross products of succeeding points
xp1=xprod(v1,v2);
xp2=xprod(v2,v3);
xp3=xprod(v3,v4);
xp4=xprod(v4,v1);
%store the cross products in array for checking
xp=[xp1,xp2,xp3,xp4];

%show the complex points generated
disp('Four random generated points:');
fprintf('Q1=%4.2f+%4.2fj\n', real(q1),imag(q1));
fprintf('Q2=%4.2f+%4.2fj\n', real(q2),imag(q2));
fprintf('Q3=%4.2f+%4.2fj\n', real(q3),imag(q3));
fprintf('Q4=%4.2f+%4.2fj\n', real(q4),imag(q4));

%plot the points
plot(real(q1),imag(q1),'Marker','*');
plot(real(q2),imag(q2),'Marker','*');
plot(real(q3),imag(q3),'Marker','*');
plot(real(q4),imag(q4),'Marker','*');
%trace the quadrilateral
plot([real(q1),real(q2)],[imag(q1),imag(q2)]);
plot([real(q2),real(q3)],[imag(q2),imag(q3)]);
plot([real(q3),real(q4)],[imag(q3),imag(q4)]);
plot([real(q4),real(q1)],[imag(q4),imag(q1)]);


%title of graph after identifying type and area
if (all(xp > 0) || all(xp <0))
    title(sprintf('CONVEX QUADRILATERAL: AREA =%4.2f',Area));%graph title
    %fprintf('CONVEX: AREA =%4.2f',Area);
elseif( any(xp ==0))
    title(sprintf('TRIANGLE: AREA =%4.2f',Area));
    %fprintf('TRIANGLE: AREA =%4.2f',Area);
else
    title(sprintf('CONCAVE QUADRILATERAL: AREA =%4.2f',Area));
    %fprintf('CONCAVE: AREA =%4.2f',Area);
end