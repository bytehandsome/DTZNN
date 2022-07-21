function [Rot, pos] = FwdKin(theta)

global a3 a4 a6 d1 d3 d5 d7;

c1 = cos(theta(1));
c2 = cos(theta(2));
c3 = cos(theta(3));
c4 = cos(theta(4));
c5 = cos(theta(5));
c6 = cos(theta(6));
c7 = cos(theta(7));

s1 = sin(theta(1));
s2 = sin(theta(2));
s3 = sin(theta(3));
s4 = sin(theta(4));
s5 = sin(theta(5));
s6 = sin(theta(6));
s7 = sin(theta(7));

%% Standard DH convention
A1 = [c1, 0, -s1, 0;
    s1, 0, c1, 0;
    0, -1, 0, d1;
    0, 0, 0, 1];

A2 = [c2, 0, s2, 0;
    s2, 0, -c2, 0;
    0, 1, 0, 0;
    0, 0, 0, 1];

A3 = [c3, 0, s3, a3*c3;
    s3, 0, -c3, a3*s3;
    0, 1, 0, d3;
    0, 0, 0, 1];

A4 = [c4, 0, -s4, a4*c4;
    s4, 0, c4, a4*s4;
    0, -1, 0, 0;
    0, 0, 0, 1];

A5 = [c5, 0, s5, 0;
    s5, 0, -c5, 0;
    0, 1, 0, d5;
    0, 0, 0, 1];

A6 = [c6, 0, -s6, a6*c6;
    s6, 0, c6, a6*s6;
    0, -1, 0, 0;
    0, 0, 0, 1];

A7 = [c7, -s7, 0, 0;...
    s7, c7, 0, 0;...
    0, 0, 1, d7;...
    0, 0, 0, 1];

H = A1*A2*A3*A4*A5*A6*A7;

pos = H(1:3, 4);
Rot = H(1:3, 1:3);

% pos(1) = d7*(c6*(s4*(s1*s3 - c1*c2*c3) + c1*c4*s2) + s6*(c5*(c4*(s1*s3 - c1*c2*c3) - c1*s2*s4) + s5*(c3*s1 + c1*c2*s3))) + d5*(s4*(s1*s3 - c1*c2*c3) + c1*c4*s2) - a4*c4*(s1*s3 - c1*c2*c3) + d3*c1*s2 - a3*s1*s3 + a6*s6*(s4*(s1*s3 - c1*c2*c3) + c1*c4*s2) - a6*c6*(c5*(c4*(s1*s3 - c1*c2*c3) - c1*s2*s4) + s5*(c3*s1 + c1*c2*s3)) + a3*c1*c2*c3 + a4*c1*s2*s4;
% pos(2) = a4*c4*(c1*s3 + c2*c3*s1) - d5*(s4*(c1*s3 + c2*c3*s1) - c4*s1*s2) - d7*(c6*(s4*(c1*s3 + c2*c3*s1) - c4*s1*s2) + s6*(c5*(c4*(c1*s3 + c2*c3*s1) + s1*s2*s4) + s5*(c1*c3 - c2*s1*s3))) + a6*c6*(c5*(c4*(c1*s3 + c2*c3*s1) + s1*s2*s4) + s5*(c1*c3 - c2*s1*s3)) + a3*c1*s3 + d3*s1*s2 - a6*s6*(s4*(c1*s3 + c2*c3*s1) - c4*s1*s2) + a3*c2*c3*s1 + a4*s1*s2*s4;
% pos(3) = d1 + d5*(c2*c4 + c3*s2*s4) - d7*(s6*(c5*(c2*s4 - c3*c4*s2) + s2*s3*s5) - c6*(c2*c4 + c3*s2*s4)) + d3*c2 + a6*s6*(c2*c4 + c3*s2*s4) - a3*c3*s2 + a4*c2*s4 + a6*c6*(c5*(c2*s4 - c3*c4*s2) + s2*s3*s5) - a4*c3*c4*s2;
% pos = transpose(pos);
%
% Rot(1,1) = c7*(s6*(s4*(s1*s3 - c1*c2*c3) + c1*c4*s2) - c6*(c5*(c4*(s1*s3 - c1*c2*c3) - c1*s2*s4) + s5*(c3*s1 + c1*c2*s3))) + s7*(s5*(c4*(s1*s3 - c1*c2*c3) - c1*s2*s4) - c5*(c3*s1 + c1*c2*s3));
% Rot(1,2) = c7*(s5*(c4*(s1*s3 - c1*c2*c3) - c1*s2*s4) - c5*(c3*s1 + c1*c2*s3)) - s7*(s6*(s4*(s1*s3 - c1*c2*c3) + c1*c4*s2) - c6*(c5*(c4*(s1*s3 - c1*c2*c3) - c1*s2*s4) + s5*(c3*s1 + c1*c2*s3)));
% Rot(1,3) = c6*(s4*(s1*s3 - c1*c2*c3) + c1*c4*s2) + s6*(c5*(c4*(s1*s3 - c1*c2*c3) - c1*s2*s4) + s5*(c3*s1 + c1*c2*s3));
% 
% Rot(2,1) = - c7*(s6*(s4*(c1*s3 + c2*c3*s1) - c4*s1*s2) - c6*(c5*(c4*(c1*s3 + c2*c3*s1) + s1*s2*s4) + s5*(c1*c3 - c2*s1*s3))) - s7*(s5*(c4*(c1*s3 + c2*c3*s1) + s1*s2*s4) - c5*(c1*c3 - c2*s1*s3));
% Rot(2,2) = s7*(s6*(s4*(c1*s3 + c2*c3*s1) - c4*s1*s2) - c6*(c5*(c4*(c1*s3 + c2*c3*s1) + s1*s2*s4) + s5*(c1*c3 - c2*s1*s3))) - c7*(s5*(c4*(c1*s3 + c2*c3*s1) + s1*s2*s4) - c5*(c1*c3 - c2*s1*s3));
% Rot(2,3) = - c6*(s4*(c1*s3 + c2*c3*s1) - c4*s1*s2) - s6*(c5*(c4*(c1*s3 + c2*c3*s1) + s1*s2*s4) + s5*(c1*c3 - c2*s1*s3));
% 
% Rot(3,1) = c7*(c6*(c5*(c2*s4 - c3*c4*s2) + s2*s3*s5) + s6*(c2*c4 + c3*s2*s4)) - s7*(s5*(c2*s4 - c3*c4*s2) - c5*s2*s3);
% Rot(3,2) = - c7*(s5*(c2*s4 - c3*c4*s2) - c5*s2*s3) - s7*(c6*(c5*(c2*s4 - c3*c4*s2) + s2*s3*s5) + s6*(c2*c4 + c3*s2*s4));
% Rot(3,3) = c6*(c2*c4 + c3*s2*s4) - s6*(c5*(c2*s4 - c3*c4*s2) + s2*s3*s5);

end