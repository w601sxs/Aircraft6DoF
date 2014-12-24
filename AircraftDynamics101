%% Aircraft Dynamics 101
clear all; close all; clc
syms x y z phi theta psi m g real 
syms Fx Fy Fz Mx My Mz real 
syms xd yd zd wx wy wz real
syms xdd ydd zdd real
syms R3 R2 R1 real
syms Ix Iy Iz real
%% Rotation

R3 = [cos(psi) sin(psi) 0;
      -sin(psi) cos(psi) 0;
      0 0 1];
R2 = [cos(theta) 0 -sin(theta);
     0 1 0;
     sin(theta) 0 cos(theta)];
R1 = [1 0 0;
      0 cos(phi) sin(phi);
      0 -sin(phi) cos(phi)];

 R=R1*R2*R3; %To rotate from body to inertial
 Rinv = R3*R2*R1;

 InertPos = Rinv*[x y z]';
 
weight = R*[0 0 -m*g]';
     
orientdot= [wy*sin(phi)*sec(theta) + wz*cos(phi)*sec(theta);
                      wy*cos(phi) - wz*sin(phi);
                      wx + wy*sin(phi)*tan(theta) + wz*cos(phi)*tan(theta)];
psid=orientdot(1); thetad=orientdot(2); phid=orientdot(3);      

%% Euler EOM
% Mx=10; My=5; Mz=3;

wxd = (Mx - (Iz - Iy)*wy*wz)/Ix;
wyd = (My - (Ix - Iz)*wx*wz)/Iy;
wzd = (Mz - (Iy - Ix)*wy*wx)/Iz;
%----------------------------------------------------------------------

%% Translation

R = [x; y; z];
w = [wx;wy;wz]; %[phid - psid*sin(theta); thetad*cos(phi) + psid*cos(theta)*sin(phi); psid*cos(theta)*cos(phi) - thetad*sin(phi)];
v = [xd; yd; zd] + cross(w,R);

a = diff(v,x)*xd + diff(v,y)*yd + diff(v,z)*zd + diff(v,xd)*xdd + diff(v,yd)*ydd + diff(v,zd)*zdd + diff(v,wx)*wxd + diff(v,wy)*wyd + diff(v,wz)*wzd;
% a = a + cross(w,v);

Rdd = solve(Fx - m*a(1), Fy - m*a(2), Fz - m*a(3),xdd, ydd, zdd);

%----------------------------------------------------------------------
S = [x y z xd yd zd phi theta psi wx wy wz]'; simplify(S);
syms vx vy vz real
xd = vx; yd=vy; zd=vz;
Rdd.xdd = subs(Rdd.xdd);
Rdd.ydd = subs(Rdd.ydd);
Rdd.zdd = subs(Rdd.zdd);
% 




