%% forward kinematics using hte screw displacement method
clc; clear; 
FontScale = 12;
LineWidth = 1.5;
LineWidth2 = 2.5;
%%%%%%%%%%
L1 = 200;
L2 = 154.75;
L3 = 318.50;
L4 = L3;
%%%%%%%%%%
theta1 = pi/6;
theta2 = pi/6;
theta3 = pi/6;
theta4 = pi/6;
theta5 = pi/6;
theta6 = pi/6;
%%%%%%%%%%
p_0 = [0; L2; -L1 - L3 - L4; 1];
% compute forward kinematics
P = forward_solve(theta1, theta2, theta3, theta4, theta5, theta6, L1, L2, L3, L4);%P = simplify(P)
%%%%%%%%%%%%%%%
f1 = figure;
hold on;
grid on;
view(3);
%%%%%%%%%%%%%%%
%frame zero coordinates
quiver3(0, 0, 0, 0, 0, L1/2, 'y', 'LineWidth',LineWidth); hold on; text(0, 0, L1/2, 'z0');
quiver3(0, 0, 0, 0, L1/2, 0, 'y', 'LineWidth',LineWidth); text(0, L1/2, 0, 'y0');
quiver3(0, 0, 0, L1/2, 0, 0, 'y', 'LineWidth',LineWidth); text(L1/2, 0, 0, 'x0');

% %initial position:
quiver3(0, 0, 0, 0, 0, -L1, 'y', 'LineWidth',LineWidth); text (0, 0, -L1/2, 'L1');
quiver3(0, 0, -L1, 0, L2, 0, 'y', 'LineWidth',LineWidth); text (0, L2/2, -L1, 'L2');

%drawing final position
quiver3(0, 0, 0, P(1), P(2), P(3), 'b', 'LineWidth',LineWidth); text (P(1)/2, P(2)/2, P(3)/2, 'position');
R1 = rotz(theta1*(180/pi)); 
R2 = roty(theta2*(180/pi)); 
R3 = rotx(theta3*(180/pi)); 
R4 = roty(theta4*(180/pi)); 
R5 = roty(theta5*(180/pi)); 
R6 = rotx(theta6*(180/pi));
%%%%%%%%%%
A = [0; 0; -L3];
p = R1*R2*R3*A; %3DOF
L3_p = quiver3(0, L2, -L1, p(1), p(2), p(3), 'r', 'LineWidth',LineWidth2);
%%%%%%%%%%
p1 = R1*R2*R3*R4*A;
L_4 = quiver3(p(1), p(2) + L2, p(3) - L1, p1(1), p1(2), p1(3), 'g', 'LineWidth',LineWidth2); %R5 and R6 are not included as they do not affect endeffector position

pu = R1*R2*R3*R4*R5*R6*[0;0;-L3];
pv = R1*R2*R3*R4*R5*R6*[0;L3;0];
pw = R1*R2*R3*R4*R5*R6*[L3;0;0];

quiver3(P(1),P(2),P(3),pu(1), pu(2), pu(3), 'y', 'LineWidth',LineWidth); text(P(1)+pu(1), P(2)+pu(2), P(3)+pu(3), 'X6');
quiver3(P(1),P(2),P(3),pv(1) , pv(2), pv(3), 'y', 'LineWidth',LineWidth); text(P(1)+pv(1), P(2)+pv(2), P(3)+pv(3), 'Y6');
quiver3(P(1),P(2),P(3),pw(1), pw(2), pw(3), 'y', 'LineWidth',LineWidth); text(P(1)+pw(1), P(2)+pw(2), P(3)+pw(3), 'Z6');

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
