function Jacob = Jacobian(theta)

global a3 a4 a6 d1 d3 d5 d7;

c1 = cos(theta(1));
c2 = cos(theta(2));
c3 = cos(theta(3));
c4 = cos(theta(4));
c5 = cos(theta(5));
c6 = cos(theta(6));
% c7 = cos(theta(7));

s1 = sin(theta(1));
s2 = sin(theta(2));
s3 = sin(theta(3));
s4 = sin(theta(4));
s5 = sin(theta(5));
s6 = sin(theta(6));
% s7 = sin(theta(7));

%% Standard DH convention
Jacob(1, 1) = d5*c1*s3*s4 - d3*s1*s2 - a3*c2*c3*s1 - a4*c1*c4*s3 - d5*c4*s1*s2 - a3*c1*s3 - a4*s1*s2*s4 + d5*c2*c3*s1*s4 - d7*c4*c6*s1*s2 + d7*c1*c6*s3*s4 + d7*c1*c3*s5*s6 - a6*c4*s1*s2*s6 + a6*c1*s3*s4*s6 - a4*c2*c3*c4*s1 - a6*c1*c3*c6*s5 - a6*c1*c4*c5*c6*s3 + d7*c2*c3*c6*s1*s4 + d7*c1*c4*c5*s3*s6 + a6*c2*c3*s1*s4*s6 + a6*c2*c6*s1*s3*s5 - a6*c5*c6*s1*s2*s4 - d7*c2*s1*s3*s5*s6 + d7*c5*s1*s2*s4*s6 - a6*c2*c3*c4*c5*c6*s1 + d7*c2*c3*c4*c5*s1*s6;
Jacob(1, 2) = c1*(d3*c2 + d5*c2*c4 - a3*c3*s2 + a4*c2*s4 + d7*c2*c4*c6 - a4*c3*c4*s2 + a6*c2*c4*s6 + d5*c3*s2*s4 + d7*c3*c6*s2*s4 - d7*c2*c5*s4*s6 + a6*c3*s2*s4*s6 + a6*c6*s2*s3*s5 - d7*s2*s3*s5*s6 + a6*c2*c5*c6*s4 - a6*c3*c4*c5*c6*s2 + d7*c3*c4*c5*s2*s6);
Jacob(1, 3) = d5*c3*s1*s4 - a3*c1*c2*s3 - a4*c3*c4*s1 - a3*c3*s1 + d5*c1*c2*s3*s4 + d7*c3*c6*s1*s4 + a6*c3*s1*s4*s6 + a6*c6*s1*s3*s5 - d7*s1*s3*s5*s6 - a4*c1*c2*c4*s3 - a6*c1*c2*c3*c6*s5 - a6*c3*c4*c5*c6*s1 + d7*c1*c2*c6*s3*s4 + d7*c1*c2*c3*s5*s6 + d7*c3*c4*c5*s1*s6 + a6*c1*c2*s3*s4*s6 - a6*c1*c2*c4*c5*c6*s3 + d7*c1*c2*c4*c5*s3*s6;
Jacob(1, 4) = a4*c1*c4*s2 - d5*c1*s2*s4 + d5*c4*s1*s3 + a4*s1*s3*s4 - d7*c1*c6*s2*s4 + d7*c4*c6*s1*s3 - a6*c1*s2*s4*s6 + a6*c4*s1*s3*s6 - d5*c1*c2*c3*c4 - a4*c1*c2*c3*s4 - d7*c1*c2*c3*c4*c6 - a6*c1*c2*c3*c4*s6 + a6*c1*c4*c5*c6*s2 - d7*c1*c4*c5*s2*s6 + a6*c5*c6*s1*s3*s4 - d7*c5*s1*s3*s4*s6 - a6*c1*c2*c3*c5*c6*s4 + d7*c1*c2*c3*c5*s4*s6;
Jacob(1, 5) = -(a6*c6 - d7*s6)*(c3*c5*s1 + c1*c2*c5*s3 + c1*s2*s4*s5 - c4*s1*s3*s5 + c1*c2*c3*c4*s5);
Jacob(1, 6) = d7*c3*c6*s1*s5 - d7*c1*c4*s2*s6 + a6*c6*s1*s3*s4 + a6*c3*s1*s5*s6 - d7*s1*s3*s4*s6 + a6*c1*c4*c6*s2 - a6*c1*c2*c3*c6*s4 + d7*c1*c2*c3*s4*s6 + d7*c1*c2*c6*s3*s5 - d7*c1*c5*c6*s2*s4 + d7*c4*c5*c6*s1*s3 + a6*c1*c2*s3*s5*s6 - a6*c1*c5*s2*s4*s6 + a6*c4*c5*s1*s3*s6 - d7*c1*c2*c3*c4*c5*c6 - a6*c1*c2*c3*c4*c5*s6;
Jacob(1, 7) = 0.0;

Jacob(2, 1) = d3*c1*s2 - a3*s1*s3 + a3*c1*c2*c3 + d5*c1*c4*s2 + a4*c1*s2*s4 - a4*c4*s1*s3 + d5*s1*s3*s4 - d5*c1*c2*c3*s4 + d7*c1*c4*c6*s2 + a6*c1*c4*s2*s6 - a6*c3*c6*s1*s5 + d7*c6*s1*s3*s4 + d7*c3*s1*s5*s6 + a6*s1*s3*s4*s6 + a4*c1*c2*c3*c4 - d7*c1*c2*c3*c6*s4 - a6*c1*c2*c3*s4*s6 - a6*c1*c2*c6*s3*s5 + a6*c1*c5*c6*s2*s4 - a6*c4*c5*c6*s1*s3 + d7*c1*c2*s3*s5*s6 - d7*c1*c5*s2*s4*s6 + d7*c4*c5*s1*s3*s6 + a6*c1*c2*c3*c4*c5*c6 - d7*c1*c2*c3*c4*c5*s6;
Jacob(2, 2) = s1*(d3*c2 + d5*c2*c4 - a3*c3*s2 + a4*c2*s4 + d7*c2*c4*c6 - a4*c3*c4*s2 + a6*c2*c4*s6 + d5*c3*s2*s4 + d7*c3*c6*s2*s4 - d7*c2*c5*s4*s6 + a6*c3*s2*s4*s6 + a6*c6*s2*s3*s5 - d7*s2*s3*s5*s6 + a6*c2*c5*c6*s4 - a6*c3*c4*c5*c6*s2 + d7*c3*c4*c5*s2*s6);
Jacob(2, 3) = a3*c1*c3 + a4*c1*c3*c4 - d5*c1*c3*s4 - a3*c2*s1*s3 - d7*c1*c3*c6*s4 - a4*c2*c4*s1*s3 - a6*c1*c3*s4*s6 - a6*c1*c6*s3*s5 + d5*c2*s1*s3*s4 + d7*c1*s3*s5*s6 + a6*c1*c3*c4*c5*c6 - d7*c1*c3*c4*c5*s6 - a6*c2*c3*c6*s1*s5 + d7*c2*c6*s1*s3*s4 + d7*c2*c3*s1*s5*s6 + a6*c2*s1*s3*s4*s6 - a6*c2*c4*c5*c6*s1*s3 + d7*c2*c4*c5*s1*s3*s6;
Jacob(2, 4) = a4*c4*s1*s2 - d5*c1*c4*s3 - a4*c1*s3*s4 - d5*s1*s2*s4 - d5*c2*c3*c4*s1 - d7*c1*c4*c6*s3 - a4*c2*c3*s1*s4 - a6*c1*c4*s3*s6 - d7*c6*s1*s2*s4 - a6*s1*s2*s4*s6 - d7*c2*c3*c4*c6*s1 - a6*c2*c3*c4*s1*s6 + a6*c4*c5*c6*s1*s2 - a6*c1*c5*c6*s3*s4 - d7*c4*c5*s1*s2*s6 + d7*c1*c5*s3*s4*s6 - a6*c2*c3*c5*c6*s1*s4 + d7*c2*c3*c5*s1*s4*s6;
Jacob(2, 5) = -(a6*c6 - d7*s6)*(c2*c5*s1*s3 - c1*c3*c5 + c1*c4*s3*s5 + s1*s2*s4*s5 + c2*c3*c4*s1*s5);
Jacob(2, 6) = a6*c4*c6*s1*s2 - d7*c1*c3*c6*s5 - a6*c1*c6*s3*s4 - a6*c1*c3*s5*s6 - d7*c4*s1*s2*s6 + d7*c1*s3*s4*s6 - d7*c1*c4*c5*c6*s3 - a6*c2*c3*c6*s1*s4 - a6*c1*c4*c5*s3*s6 + d7*c2*c3*s1*s4*s6 + d7*c2*c6*s1*s3*s5 - d7*c5*c6*s1*s2*s4 + a6*c2*s1*s3*s5*s6 - a6*c5*s1*s2*s4*s6 - d7*c2*c3*c4*c5*c6*s1 - a6*c2*c3*c4*c5*s1*s6;
Jacob(2, 7) = 0.0;

Jacob(3, 1) = 0.0;
Jacob(3, 2) = d5*c2*c3*s4 - a3*c2*c3 - d5*c4*s2 - a4*s2*s4 - a4*c2*c3*c4 - d3*s2 - d7*c4*c6*s2 - a6*c4*s2*s6 + d7*c2*c3*c6*s4 + a6*c2*c3*s4*s6 + a6*c2*c6*s3*s5 - a6*c5*c6*s2*s4 - d7*c2*s3*s5*s6 + d7*c5*s2*s4*s6 - a6*c2*c3*c4*c5*c6 + d7*c2*c3*c4*c5*s6;
Jacob(3, 3) = -s2*(d5*s3*s4 - a4*c4*s3 - a3*s3 - a6*c3*c6*s5 + d7*c6*s3*s4 + d7*c3*s5*s6 + a6*s3*s4*s6 + d7*c4*c5*s3*s6 - a6*c4*c5*c6*s3);
Jacob(3, 4) = a4*c2*c4 - d5*c2*s4 + d5*c3*c4*s2 - d7*c2*c6*s4 + a4*c3*s2*s4 - a6*c2*s4*s6 + d7*c3*c4*c6*s2 - d7*c2*c4*c5*s6 + a6*c3*c4*s2*s6 + a6*c2*c4*c5*c6 + a6*c3*c5*c6*s2*s4 - d7*c3*c5*s2*s4*s6;
Jacob(3, 5) = (a6*c6 - d7*s6)*(c5*s2*s3 - c2*s4*s5 + c3*c4*s2*s5);
Jacob(3, 6) = a6*c2*c4*c6 - d7*c2*c4*s6 - d7*c2*c5*c6*s4 + a6*c3*c6*s2*s4 - a6*c2*c5*s4*s6 - d7*c3*s2*s4*s6 - d7*c6*s2*s3*s5 - a6*s2*s3*s5*s6 + d7*c3*c4*c5*c6*s2 + a6*c3*c4*c5*s2*s6;
Jacob(3, 7) = 0.0;

Jacob(4, 1) = 0.0;
Jacob(4, 2) = -s1;
Jacob(4, 3) = c1*s2;
Jacob(4, 4) = c3*s1 + c1*c2*s3; 
Jacob(4, 5) = s4*(s1*s3 - c1*c2*c3) + c1*c4*s2; 
Jacob(4, 6) = c5*(c3*s1 + c1*c2*s3) - s5*(c4*(s1*s3 - c1*c2*c3) - c1*s2*s4); 
Jacob(4, 7) = s6*(s5*(c3*s1 + c1*c2*s3) + c5*(c4*(s1*s3 - c1*c2*c3) - c1*s2*s4)) + c6*(s4*(s1*s3 - c1*c2*c3) + c1*c4*s2); 

Jacob(5, 1) = 0.0;
Jacob(5, 2) = c1;
Jacob(5, 3) = s1*s2;
Jacob(5, 4) = c2*s1*s3 - c1*c3; 
Jacob(5, 5) = c4*s1*s2 - s4*(c1*s3 + c2*c3*s1); 
Jacob(5, 6) = s5*(c4*(c1*s3 + c2*c3*s1) + s1*s2*s4) - c5*(c1*c3 - c2*s1*s3); 
Jacob(5, 7) = - s6*(s5*(c1*c3 - c2*s1*s3) + c5*(c4*(c1*s3 + c2*c3*s1) + s1*s2*s4)) - c6*(s4*(c1*s3 + c2*c3*s1) - c4*s1*s2); 

Jacob(6, 1) = 1.0;
Jacob(6, 2) = 0.0;
Jacob(6, 3) = c2;
Jacob(6, 4) = -s2*s3;
Jacob(6, 5) = c2*c4 + c3*s2*s4; 
Jacob(6, 6) = s5*(c2*s4 - c3*c4*s2) - c5*s2*s3; 
Jacob(6, 7) = c6*(c2*c4 + c3*s2*s4) - s6*(c5*(c2*s4 - c3*c4*s2) + s2*s3*s5);

end