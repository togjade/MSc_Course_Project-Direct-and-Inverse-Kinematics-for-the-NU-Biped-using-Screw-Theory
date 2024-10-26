%% correct
% compute forward kinematics to substitute their values into the inverse
% kinematics
clc; clear; 
syms theta1 theta2 theta3 theta4 theta5 theta6 L1 L2 L3 L4 p_x p_y p_z w_x w_y w_z v_x v_y v_z u_x u_y u_z
L1 = 200;
L2 = 154.75;
L3 = 318.50;
L4 = L3;
%%%%%%%%%%
theta1 = 0; theta2 = 0;
theta3 = pi/6; theta4 = 0;
theta5 = 0; theta6 = 0;
%%%%%%%%%%
R1 = rotz(theta1*(180/pi)); 
R2 = roty(theta2*(180/pi)); 
R3 = rotx(theta3*(180/pi)); 
R4 = roty(theta4*(180/pi)); 
R5 = roty(theta5*(180/pi)); 
R6 = rotx(theta6*(180/pi));
p_0 = [0; L2; -L1 - L3 - L4; 1];
% compute forward kinematics
P = forward_solve(theta1, theta2, theta3, theta4, theta5, theta6, L1, L2, L3, L4);%P = simplify(P)
u = R1*R2*R3*R4*R5*R6*[0;0;-1];
v = R1*R2*R3*R4*R5*R6*[0;1;0];
w = R1*R2*R3*R4*R5*R6*[1;0;0];
%%%%%%%%%%%%%%%%%%%%
p_x =  P(1,1);
p_y = P(2,1);
p_z =  P(3,1);

w_x = w(1,1);
w_y = w(2,1);
w_z = w(3,1);

v_x = v(1,1);
v_y = v(2,1);
v_z = v(3,1);

u_x = u(1,1);
u_y = u(2,1);
u_z = u(3,1);

NewP = NewPoint(p_x, p_y, p_z, u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z);

u_x = NewP(1,1);
u_y=  NewP(2,1);
u_z =  NewP(3,1);

v_x = NewP(1,2);
v_y =  NewP(2,2);
v_z =  NewP(3,2);

w_x = NewP(1,3);
w_y =  NewP(2,3);
w_z =  NewP(3,3);

p_x = NewP(1,4); 
p_y =  NewP(2,4);
p_z =  NewP(3,4);

thetas = IKSol(p_x, p_y, p_z, v_x, v_y, v_z, w_x, w_y, w_z);

function C = NewPoint(p_x, p_y, p_z, u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z)
L1= 200;
L2 = 154.75;

A = [u_x v_x w_x p_x;
     u_y v_y w_y p_y;
     u_z v_z w_z p_z;
     0    0   0  1];

T1 = [1 0 0 0;
      0 1 0 +L2;
      0 0 1 -L1;
      0 0 0 1];

C = inv(A)*T1;

end
function B = IKSol(p_x, p_y, p_z, v_x, v_y, v_z, w_x, w_y, w_z)
L3 =318.5;
L4=L3;
w_0 = [-1; 0; 0];
v_0 = [0; 1; 0];
v = [v_x; v_y; v_z];
w = [w_x; w_y; w_z];
P_0 = [- L3 - L4; 0 ; 0; 1];
P = [p_x; p_y; p_z; 1];


% theta6
theta6 = atan(p_y/p_x);
theta6d = pi+theta6;

k1 = p_x*cos(theta6)+p_y*sin(theta6);

%theta4
theta4 = acos((p_z^2+k1^2-L3^2-L4^2)/(2*L3*L4));
theta4d = -theta4;

%theta5

sin5 = (L3*(p_z*cos(theta4) + k1*sin(theta4)) + p_z*L4)/(L3^2 + L4^2 + 2*L3*L4*cos(theta4));
cos5 = (-k1 + L3*sin(theta4)*sin5)/(L3*cos(theta4)+L4);
theta5 = atan2(sin5, cos5);

%theta3
k3 = w_y*cos(theta6)-w_x*sin(theta6);
k2 = w_x*cos(theta4+theta5)*cos(theta6) - w_z*sin(theta4+theta5) + w_y*cos(theta4+theta5)*sin(theta6);
theta3 = atan(k3/k2);


%theta2
k4 = w_z*cos(theta4+theta5) + w_x*sin(theta4+theta5)*cos(theta6) + w_y*sin(theta4+theta5)*sin(theta6);
theta2 = atan2(k4, -k2/cos(theta3));


%theta1
k5 = v_x*cos(theta4+theta5)*cos(theta6) - v_z*sin(theta4+theta5) + v_y*cos(theta4+theta5)*sin(theta6);

k7 = v_z*cos(theta4+theta5) + v_x*sin(theta4+theta5)*cos(theta6) + v_y*sin(theta4+theta5)*sin(theta6);

sin1 = -k7/cos(theta2);

cos1 = (k7*k4*cos(theta3)-k5*cos(theta2))/sin(theta3);
theta1 = atan2(sin1, cos1);
%theta1d = atan2(sin1, -sqrt(1-sin1^2));

B = [theta1 theta2 theta3 theta4 theta5 theta6];

end 
%%
function A = make_matrix(s_x, s_y, s_z, theta, s_ox, s_oy, s_oz)
a_11 = (s_x^2 - 1)*(1 - cos(theta)) + 1;
a_12 = s_x*s_y*(1 - cos(theta)) - s_z*sin(theta);
a_13 = s_x*s_z*(1 - cos(theta)) + s_y*sin(theta);
a_21 = s_y*s_x*(1 - cos(theta)) + s_z*sin(theta);
a_22 = (s_y^2 - 1)*(1 - cos(theta)) + 1;
a_23 = s_y*s_z*(1 - cos(theta)) - s_x*sin(theta);
a_31 = s_z*s_x*(1 - cos(theta)) - s_y*sin(theta);
a_32 = s_z*s_y*(1 - cos(theta)) + s_x*sin(theta);
a_33 = (s_z^2 - 1)*(1 - cos(theta)) + 1;
a_14 = -s_ox*(a_11 - 1) - s_oy*a_12 - s_oz*a_13;
a_24 = -s_ox*a_21 - s_oy*(a_22 - 1) - s_oz*a_23;
a_34 = -s_ox*a_31 - s_oy*a_32 - s_oz*(a_33 - 1);

A = [a_11 a_12 a_13 a_14;
     a_21 a_22 a_23 a_24;
     a_31 a_32 a_33 a_34;
     0    0    0    1];
 
end 


function P = forward_solve(theta1, theta2, theta3, theta4, theta5, theta6, L1, L2, L3, L4)
P_0 = [0; L2; -L1 - L3 - L4; 1];

A_01 = make_matrix(0, 0, 1, theta1, 0, L2, -L1);
A_12 = make_matrix(0, 1, 0, theta2, 0, L2, -L1);
A_23 = make_matrix(1, 0, 0, theta3, 0, L2, -L1);

A_34 = make_matrix(0, 1, 0, theta4, 0, L2, -L1 - L3);

A_45 = make_matrix(0, 1, 0, theta5, 0, L2, -L1 - L3 - L4);
A_56 = make_matrix(1, 0, 0, theta6, 0, L2, -L1 - L3 - L4);


A_06 = A_01*A_12*A_23*A_34*A_45*A_56;

P = A_06*P_0;
end 
