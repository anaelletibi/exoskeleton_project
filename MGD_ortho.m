function out1 = MGD_ortho(q)
% MGD du robot à axes orthogonaux : T0j, j=1..7
% usage : T01 = out1(:,:,1) à T07 = out1(:,:,7)
% T0j matrice homogène 4x4


q1=q(1);
q2=q(2);
q3=q(3);
q4=q(4);
q5=q(5);
q6=q(6);
q7=q(7);


t1=q1;
ct1=cos(t1);
st1=sin(t1);
d1=0;
a1= pi/2;
r1=0;
T01=[ct1 -st1 0 d1; ...
cos(a1)*st1 cos(a1)*ct1 -sin(a1) -r1*sin(a1); ...
sin(a1)*st1 sin(a1)*ct1 cos(a1) r1*cos(a1); ...
0 0 0 1];

t2=q2-pi/2;
ct2=cos(t2);
st2=sin(t2);
d2=0;
a2= -pi/2;
ca2= cos(a2);
sa2=sin(a2);
r2=0;
T12=[ct2 -st2 0 d2; ...
ca2*st2 ca2*ct2 -sa2 -r2*sa2; ...
sa2*st2 sa2*ct2 ca2 r2*ca2; ...
0 0 0 1];
T02=T01*T12;

t3=q3+pi/2;
ct3=cos(t3);
st3=sin(t3);
d3=0;
a3= -pi/2;
ca3= cos(a3);
sa3=sin(a3);
r3=0;
T23=[ct3 -st3 0 d3; ...
ca3*st3 ca3*ct3 -sa3 -r3*sa3; ...
sa3*st3 sa3*ct3 ca3 r3*ca3; ...
0 0 0 1];
T03=T02*T23;

t4=q4-pi;
ct4=cos(t4);
st4=sin(t4);
d4=0.300;
a4= 0;
ca4= cos(a4);
sa4=sin(a4);
r4=0;
T34=[ct4 -st4 0 d4; ...
ca4*st4 ca4*ct4 -sa4 -r4*sa4; ...
sa4*st4 sa4*ct4 ca4 r4*ca4; ...
0 0 0 1];
T04=T03*T34;

t5=q5;
ct5=cos(t5);
st5=sin(t5);
d5=0;
a5= -pi/2;
ca5= cos(a5);
sa5=sin(a5);
r5=0.230;
T45=[ct5 -st5 0 d5; ...
ca5*st5 ca5*ct5 -sa5 -r5*sa5; ...
sa5*st5 sa5*ct5 ca5 r5*ca5; ...
0 0 0 1];
T05=T04*T45;

t6=q6+pi/2;
ct6=cos(t6);
st6=sin(t6);
d6=0;
a6= pi/2;
ca6= cos(a6);
sa6=sin(a6);
r6=0;
T56=[ct6 -st6 0 d6; ...
ca6*st6 ca6*ct6 -sa6 -r6*sa6; ...
sa6*st6 sa6*ct6 ca6 r6*ca6; ...
0 0 0 1];
T06=T05*T56;

t7=q7;
ct7=cos(t7);
st7=sin(t7);
d7=0;
a7= pi/2;
ca7= cos(a7);
sa7=sin(a7);
r7=0;
T67=[ct7 -st7 0 d7; ...
ca7*st7 ca7*ct7 -sa7 -r7*sa7; ...
sa7*st7 sa7*ct7 ca7 r7*ca7; ...
0 0 0 1];
T07=T06*T67;

out1(:,:,1) = T01;
out1(:,:,2) = T02;
out1(:,:,3) = T03;
out1(:,:,4) = T04;
out1(:,:,5) = T05;
out1(:,:,6) = T06;
out1(:,:,7) = T07;
