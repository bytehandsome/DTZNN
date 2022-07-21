clc;
clearvars;
close all;

global T radius zeta;
global a3 a4 a6 d1 d3 d5 d7;

% Add paths
path1 = genpath('vrep_utils');
path2 = genpath('robot_kinematics');
path3 = genpath('plot_tools');
addpath(path1, path2, path3);
% OSQP
% Ubuntu
% addpath '/media/wellbeinglwb/Work/Research/CommonTools/osqp/osqp-0.5.0-matlab-linux64'
% Windows 10
% addpath 'G:\Research\CommonTools\osqp\osqp-0.5.0-matlab-windows64'

%% Parameters for simulation
T = 20;
tau = 0.001;
simuTime = 0:tau:T;
iterNum = length(simuTime);
radius = 0.1; % 0.3 for limit 1.5

% Robot parameters
% d1 = 0.333; % With respect to the robot base
d1 = 1.053; % With respect to the world frame if the robot is attached to a table, 0.333 + 0.72
a3 = 0.0825;
d3 = 0.316;
a4 = -0.0825;
d5 = 0.384;
a6 = -0.088;
% d7 = 0.107; % Without hand
d7 = 0.2194; % With hand, 0.107 + 0.1034 + 0.009


% Define desired position of end effector
% pOT_Od = pOT_O + dpOT_O;
pOT_Od = [ 0.5793;0;0.9043]; % Be compatibel with V-REP
% pOT_Od = [0.6; 0.2; 1.0];
% Define desired orientation of end effector
% rOTd = R2q(ROT*RT0Td);
% ROTd = [0, 1, 0;...
%     1, 0, 0;...
%     0, 0, -1];
% ROTd = [-0.7071 ,   0.7071 ,   0.0000;
%     0.7071,    0.7071   ,      0;
%    -0.0000,    0.0000 ,  -1.0000];
ROTd =  [     0.8192,    0.5736 ,   0.0000;
    0.5736  , -0.8192 ,  -0.0000;
         0  ,  0.0000 ,  -1.0000];
% ROTd = ROT;
rOTd = R2q(ROTd); % Be compatibel with V-REP


numJoints = 7;
theta = [0.0; pi/12.0; 0.0; -2.0*pi/3.0; 0.0; 3.0*pi/4.0; pi/4.0]; % Joint states in V-REP
jointOffset = [0.0; 0.0; 0.0; 0.0; 0.0; -pi; 0.0];
% q0 = theta + jointOffset; % Joint states in DH sense
q0 =  [   28.6752;16.7063;-11.3103; -102.4756; 3.6896;118.8248-180;106.0499]*pi/180;

% Joint angle limits (rad) in V-REP
qMax = [166; 101; 166; -4; 166; 215; 166]*pi/180;
qMin = [-166; -101; -166; -176; -166; -1; -166]*pi/180;

% Joint angle limits (rad) in DH sense
qMax = qMax + jointOffset;
qMin = qMin + jointOffset;

% Joint velocity limits (rad/s) in V-REP
dqMax =  [150; 150; 150; 150; 180; 180; 180]*pi/180;
dqMin = -dqMax;

%% Neural algorithm parameters
limit = 0.95;
% Proportional control gain
% ncp1 best 3 3 1000 5
Kp = 3; %5;
% Derivative control gain
Kd = 3; %5;
gamma = 1000;  %10; % Convergence parameter tau*gamma = 0.5  1000
lambda = 5; %5; % Convergence parameter of orientation  5
zeta = 0; % Control feedback
rho = 1; % Inequality scaling factor
kappa = 6; % Joint limit avoidance
% For selection of kappa, see 11j_rcm.pdf or 03j_IEEEtnn:
% 2*max(dqMax./(qMax - qMin), -dqMin./(qMax - qMin))
% max((dqMax - dqMin)./(qMax - qMin))

%% Variables
q = q0;
dq = zeros(numJoints, 1);
Jacob = Jacobian(q);
dotJacob = dotJacobian(q, dq);
Jp = Jacob(1:3, :);
% Jw = Jacob(4:6, :);
% dotA = dotJacob(1:3, :);
% dotJw = dotJacob(4:6, :);
[m, n] = size(Jp);
C = [eye(n, n); -eye(n, n)];
p = size(C, 1);

x = zeros(numJoints + m + p, iterNum + 1);
jointAngle = zeros(numJoints, iterNum + 1);
jointVelocity = zeros(numJoints, iterNum + 1);
jointAngle(:, 1) = q; % DH sense

% Initilize data for recording trace
pOT_trace = zeros(3, iterNum);
rOT_trace = zeros(4, iterNum);
Err_trace = zeros(numJoints, iterNum);
q_trace = zeros(numJoints, iterNum);

%% Update law
tic;
for i = 1:iterNum
    disp(i);
    
    [ROT, pOT_O] = FwdKin(q);
    % Convert rotation into quaternion
    rOT = R2q(ROT);
    
    pOT_trace(:, i) = pOT_O;
    rOT_trace(:, i) = rOT;
    q_trace(:, i) = q;
    
    Jacob = Jacobian(q);
    dotJacob = dotJacobian(q, dq);
    Jp = Jacob(1:3, :);
    Jw = Jacob(4:6, :);
    dotJp = dotJacob(1:3, :);
    dotJw = dotJacob(4:6, :);
    
    u = dq;
    v = x(n+1:m+n, i);
    w = x(m+n+1:m+n+p, i);
    
    A = Jp;
    dotA = dotJp;
    vecB = Kp*(pOT_Od - pOT_O);
    dotB = Kp*(pOT_Od - pOT_O) - Kd*Jp*dq;
    
    %% TSMC method
    c1 = qMax*limit;
    c2 = qMax - c1;
    c3 = qMin*limit;
    c4 = qMin - c3;
    
    dqP = dqMax - dqMax.*(sin(0.5*pi.*(sin(0.5*pi.*(q - c1)./c2)).^2)).^2;
    dqM = dqMin - dqMin.*(sin(0.5*pi.*(sin(0.5*pi.*(q - c3)./c4)).^2)).^2;
    
    ddqP = -(pi^2*cos((pi.*(c1 - q))./(2*c2)).*sin((pi.*(c1 - q))./(2*c2)).*cos((pi.*sin((pi.*(c1 - q))./(2*c2)).^2)./2).*sin((pi*sin((pi.*(c1 - q))./(2*c2)).^2)./2))./c2;
    ddqP = - dqMax.*ddqP.*dq;
    ddqM = -(pi^2*cos((pi.*(c3 - q))./(2*c4)).*sin((pi.*(c3 - q))./(2*c4)).*cos((pi.*sin((pi.*(c3 - q))./(2*c4)).^2)./2).*sin((pi*sin((pi.*(c3 - q))./(2*c4)).^2)./2))./c4;
    ddqM = - dqMin.*ddqM.*dq;
    
    II = find(q < c1);
    dqP(II) = dqMax(II);
    ddqP(II) = 0;
    
    II = find(q > c3);
    dqM(II) = dqMin(II);
    ddqM(II) = 0;

    %% TNNLS method
    %     c1 = qMax*limit;
    %     c2 = qMax - c1;
    %     c3 = qMin*limit;
    %     c4 = qMin - c3;
    %
    %     tempExp1 = exp(-12 * (q - c3) ./ (qMin - c3) + 5);
    %     tempExp2 = exp(-12 * (q - c1) ./ (qMax - c1) + 5);
    %
    %     dqM = - dqMin ./ (1 + tempExp1) + dqMin;
    %     dqP = - dqMax ./ (1 + tempExp2) + dqMax;
    %
    %     ddqM = (12*dqMin.*tempExp1)./((tempExp1 + 1).^2.*(c3 - qMin));
    %     ddqM = ddqM.*dq;
    %
    %     ddqP = (12*dqMax.*tempExp2)./((tempExp2 + 1).^2.*(c1 - qMax));
    %     ddqP = ddqP.*dq;
    
    %% Conventional method
    % dqM = max(kappa*(qMin - q), dqMin);
    % dqP = min(kappa*(qMax - q), dqMax);
    
    vecD = [dqP; -dqM];
    dotD = [ddqP; -ddqM];
    
    % Error vector defined by Daniel Kruse
    er = quatproduct(rOTd)*quatcomplement(rOT); % rOTd x rOT*
    er(1) = 1 - er(1);
    Err_trace(:,i) = [pOT_Od - pOT_O; er];
    % Jr is the E(q) in the paper: er is used, why not rOT?
    % Jr = quatjacobian(er);
    Jr = quatjacobian(rOT); % E(q) in the paper
    Jrv = Jr*Jw; % Jq = E * Jo in the paper
    
    dotrOT = Jr*Jw*dq;
    % doter = quatproduct(rOTd)*quatcomplement(dotrOT); % dotrOTd = 0
    % dotJr = quatjacobian(doter);
    dotJr = quatjacobian(dotrOT);
    dotJrv = dotJr*Jw + Jr*dotJw;
    
    Jb = quatjacobianB(rOTd); % B in the paper
    matrixA = Jb*Jrv; % A = B * Jq in the paper
    
    %     W = eye(n, n);
    %     dotW = zeros(n, n);
    %     vecQ = zeros(n, 1);
    %     dotQ = zeros(n, 1);
    
    W = transpose(matrixA)*matrixA;
    vecQ = lambda*transpose(matrixA)*Jb*rOT;
    
    matrixAdot = Jb*dotJrv;
    dotW = transpose(matrixAdot)*matrixA + transpose(matrixA)*matrixAdot;
    dotQ = lambda*transpose(matrixAdot)*Jb*rOT + lambda*transpose(matrixA)*Jb*dotrOT; % assume dotJb = 0
    
    dotC = zeros(size(C));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Perturbed Fisher-Burmeister function: a + b - sqrt(a.*a + b.*b + delta)
%         e = vecD - C*u;
%         feps = 0.0001;
%         D1 = diag(e);
%         D2 = diag(w);
%     
%         z = D1*e + D2*w;
%         z = sqrt(z + feps);
%     
%         D3 = diag(e./z);
%         D4 = diag(w./z);
%     
%         M = [W, transpose(A), transpose(C);...
%             A, zeros(m, m), zeros(m, p);...
%             -C + D3*C, zeros(p, m), eye(size(D4))-D4];
%     
%         P = [dotW, transpose(dotA), transpose(dotC);...
%             dotA, zeros(m, m), zeros(m, p);...
%             -dotC + D3*dotC, zeros(p, m), zeros(p, p)];
%     
%         vecDotV = [dotQ; -dotB; dotD - D3*dotD];
%     
%         err = [W*u + vecQ + transpose(A)*v + transpose(C)*w; A*u - vecB; w + e - z];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% NCP function I: 0.5*min(0, a+b)*min(0, a+b) - a.*b
        e = vecD - C*u;
        eta = min(0, w + e);
        D1 = diag(w);
        D2 = diag(e);
        D3 = diag(eta);
        delta = 0.5; % ADD by whm
    
        M = [W, transpose(A), transpose(C);...
            A, zeros(m, m), zeros(m, p);...
            (D1 - D3)*C, zeros(p, m), D3 - D2];
%         M = [W, transpose(A), transpose(C);...
%             A, zeros(m, m), zeros(m, p);...
%             (D1 - D3)*C-delta*(C), zeros(p, m), (D3 - D2)+delta*eye(14)];
    
        P = [dotW, transpose(dotA), transpose(dotC);...
            dotA, zeros(m, m), zeros(m, p);...
            (D1 - D3)*dotC, zeros(p, m), zeros(p, p)];
    
        vecDotV = [dotQ; -dotB; (D3 - D1)*dotD];
    
        err = [W*u + vecQ + transpose(A)*v + transpose(C)*w; A*u - vecB; 0.5 * eta.^2 - D1*e];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% NCP function II: min(0, a)*min(0, a) + min(0, b)*min(0, b) - a.*b
    %     e = vecD - C*u;
    %     eta1 = min(0, w);
    %     eta2 = min(0, e);
    %     D1 = diag(w);
    %     D2 = diag(e);
    %     D3 = diag(eta1);
    %     D4 = diag(eta2);
    %
    %     M = [W, transpose(A), transpose(C);...
    %         A, zeros(m, m), zeros(m, p);...
    %         (D1 - 2*D4)*C, zeros(p, m), 2*D3 - D2];
    %
    %     P = [dotW, transpose(dotA), transpose(dotC);...
    %         dotA, zeros(m, m), zeros(m, p);...
    %         (D1 - 2*D4)*dotC, zeros(p, m), zeros(p, p)];
    %
    %     vecDotV = [dotQ; -dotB; (2*D4 - D1)*dotD];
    %
    %     err = [W*u + vecQ + transpose(A)*v + transpose(C)*w; A*u - vecB; eta1.^2 + eta2.^2 - D1*e];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% NCP function III (MaxQ method): max(0, -e)*max(0, -e)*max(0, -e)
%     e = vecD - C*u;
%     eta = max(0, -e);
%     D1 = diag(w);
%     D2 = diag(eta);
%     D3 = D1*D2;
%     D4 = D2*D2;
%     
%     M = [W + 6*transpose(C)*D3*C, transpose(A), 3*transpose(C)*D4;...
%         A, zeros(m, m), zeros(m, p);...
%         3*D4*C, zeros(p, m), zeros(p, p)];
%     
%     P = [dotW + 6*transpose(C)*D3*dotC, transpose(dotA), 3*transpose(dotC)*D4;...
%         dotA, zeros(m, m), zeros(m, p);...
%         3*D4*dotC, zeros(p, m), zeros(p, p)];
%     
%     vecDotV = [dotQ - 6*transpose(C)*D3*dotD; -dotB; -3*D4*dotD];
%     
%     err = [W*u + vecQ + transpose(A)*v + 3*transpose(C)*D1*D2*eta; A*u - vecB; D2*D2*eta];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    x(:, i + 1) = x(:, i) - tau*pinv(M)*(P*x(:, i) + vecDotV + gamma*err);
    
    dqTemp = x(1:numJoints, i + 1);
    q = q + tau*dq + 0.5*tau*(dqTemp - dq);
    dq = dqTemp;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% OSQP algorithm
    % min. x'Wx/2 + q'x
    % s.t. l <= Ex <= u
    
    %     E = [A; eye(7, 7)];
    %     % lower_osqp = [vecB; dqM];
    %     % upper_osqp = [vecB; dqP];
    %     lower_osqp = [vecB; -10 * ones(7, 1)];
    %     upper_osqp = [vecB; 10 * ones(7, 1)];
    %     % E = [A; Jrv(2:4,:); eye(7, 7)];
    %     % lower_osqp = [vecB; er(2:4); -10 * ones(7, 1)];
    %     % upper_osqp = [vecB; er(2:4); 10 * ones(7, 1)];
    %     sparseW = sparse(W);
    %     sparseE = sparse(E);
    %
    %     % Create an OSQP object
    %     prob = osqp;
    %
    %     % Setup workspace and change alpha parameter
    %     prob.setup(sparseW, vecQ, sparseE, lower_osqp, upper_osqp, 'eps_abs', 1.0e-16, 'eps_rel', 1.0e-16, 'max_iter', 1000,  'verbose', 1);
    %
    %     % Solve problem
    %     output = prob.solve();
    %     dq = output.x;
    %
    %     % disp(Jb *rOT - er(2:4))
    %     q = q + tau*dq;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    jointAngle(:, i + 1) = q; % DH sense
    jointVelocity(:, i + 1) = dq;
end
toc;

clear t q dq
t = simuTime;
q = jointAngle;
dq = jointVelocity;

%% VREP visualization in an off-line manner
disp('Program started');
vrep = remApi('remoteApi');
vrep.simxFinish(-1);
clientID = vrep.simxStart('127.0.0.1', 19999, true, true, 5000, 5);

if (clientID > -1)
    disp('Connected to remote API server');
    
    % Synchronous mode
    vrep.simxSynchronous(clientID, true);
    
    % Get robot handles
    for i = 1:numJoints
        motorName = ['Franka_joint' num2str(i)];
        [returnCode, handle] = vrep.simxGetObjectHandle(clientID, motorName, vrep.simx_opmode_blocking);
        if returnCode ~= vrep.simx_return_ok
            disp(['Handle of motor ' motorName ' is not found!']);
        else
            disp(['Find motor ' motorName]);
        end
    end
    
    % Start simulation
    vrep.simxStartSimulation(clientID, vrep.simx_opmode_oneshot);
    pause(3);
    
    i = 1;
    while (vrep.simxGetConnectionId(clientID) ~= -1 && i < length(t))
        theta = q(:, i) - jointOffset;
        [returnCode, ~, ~, ~, ~] = vrep.simxCallScriptFunction(clientID, 'Franka', vrep.sim_scripttype_childscript, 'setJointPos_function', [], theta, '', [], vrep.simx_opmode_blocking);
        
        if (returnCode == vrep.simx_return_ok)
            disp(['Call remote function successfully and i is ' num2str(i)]);
        else
            disp('Failed to call remote function!');
        end
        if i == 1
            % Set initial states of the simulated robot in V-REP
            disp('Set initial states of the simulated robot in V-REP');
            pause(0.001);
        else
            pause(0.001);
        end
        i = i + 10; % i = i + 1; Using 10 for quick visualization
        
        % NOTE: Put the update law (inside the previous for loop) here for 
        % on-line visualization, refer to the other programs
         
        % Trigger the next simulation step
        vrep.simxSynchronousTrigger(clientID);
    end
    vrep.simxStopSimulation(clientID, vrep.simx_opmode_oneshot);
    [returnCode, pingTime] = vrep.simxGetPingTime(clientID);
    vrep.simxFinish(clientID);
    disp('Connection is colosed!')
else
    disp('Failed connecting to remote API server');
end
vrep.delete();
disp('Task is finished and program ended');

% save data.mat simuTime jointAngle jointVelocity

%% Data processing
jointAngle = jointAngle';
jointVelocity = jointVelocity';
jointAngle = jointAngle(1:end-1, :); % discard the last one
jointVelocity = jointVelocity(1:end-1, :); % discard the last one

% [~, pos] = FwdKin(q0);
% ix = pos(1);
% iy = pos(2);
% iz = pos(3);
% t = simuTime;
% phi_sin=2*pi*sin(0.5*pi*t/T);
% phi=phi_sin.*sin(0.5*pi*t/T);
% phiDot=phi_sin*pi.*cos(0.5*pi*t/T)/T;
% 
% rx=radius*cos(phi).*(1-cos(phi))+ix;
% ry=radius*sin(phi).*(1-cos(phi))+iy;
% rz=iz*ones(size(rx));

P = zeros(3, length(t));
for i = 1:length(t)
    [~, P(:, i)] = FwdKin(jointAngle(i, :));
end

% figure
% plot(rx, ry, 'r--');
% hold on;
% plot(P(1, :), P(2, :), 'b-.');

figure(2)
subplot(1,2,1);
plot(t,Err_trace);
title('err');
% figure
% plot(t, q_trace);
% figure
% plot(t, pOT_trace);
% figure
% plot(t, rOT_trace);

% figure
% plot(t, jointAngle*180/pi); % DH sense

subplot(1,2,2);
plot(t, jointVelocity);
title('jointVelocity');

% figure
% % plot3(rx, ry, rz, 'r--');
% % hold on;
% plot3_fewer_markers(P(1, :), P(2, :), P(3, :), 20, '-bp', 'MarkerFaceColor', 'm', 'MarkerSize', 6, 'LineWidth', 1);
% zlim([0, 2])

% figure
% plot(t, rx - P(1, :), 'r');
% hold on;
% plot(t, ry - P(2, :), 'b--');
% hold on;
% plot(t, rz - P(3, :), 'c-.');

jointLimit = zeros(size(jointAngle, 1), 2*numJoints);
for index = 1:size(jointAngle, 1)
    q = jointAngle(index, :);
    q = q';
    dq = jointVelocity(index, :);
    dq = dq';
    
    %% TSMC method
    c1 = qMax*limit;
    c2 = qMax - c1;
    c3 = qMin*limit;
    c4 = qMin - c3;
    
    dqP = dqMax - dqMax.*(sin(0.5*pi.*(sin(0.5*pi.*(q - c1)./c2)).^2)).^2;
    dqM = dqMin - dqMin.*(sin(0.5*pi.*(sin(0.5*pi.*(q - c3)./c4)).^2)).^2;
    
    ddqP = -(pi^2*cos((pi.*(c1 - q))./(2*c2)).*sin((pi.*(c1 - q))./(2*c2)).*cos((pi.*sin((pi.*(c1 - q))./(2*c2)).^2)./2).*sin((pi*sin((pi.*(c1 - q))./(2*c2)).^2)./2))./c2;
    ddqP = - dqMax.*ddqP.*dq;
    ddqM = -(pi^2*cos((pi.*(c3 - q))./(2*c4)).*sin((pi.*(c3 - q))./(2*c4)).*cos((pi.*sin((pi.*(c3 - q))./(2*c4)).^2)./2).*sin((pi*sin((pi.*(c3 - q))./(2*c4)).^2)./2))./c4;
    ddqM = - dqMin.*ddqM.*dq;
    
    II = find(q < c1);
    dqP(II) = dqMax(II);
    ddqP(II) = 0;
    
    II = find(q > c3);
    dqM(II) = dqMin(II);
    ddqM(II) = 0;
    
    %% TNNLS method
    %     c1 = qMax*limit;
    %     c2 = qMax - c1;
    %     c3 = qMin*limit;
    %     c4 = qMin - c3;
    %
    %     tempExp1 = exp(-12 * (q - c3) ./ (qMin - c3) + 5);
    %     tempExp2 = exp(-12 * (q - c1) ./ (qMax - c1) + 5);
    %
    %     dqM = - dqMin ./ (1 + tempExp1) + dqMin;
    %     dqP = - dqMax ./ (1 + tempExp2) + dqMax;
    %
    %     ddqM = (12*dqMin.*tempExp1)./((tempExp1 + 1).^2.*(c3 - qMin));
    %     ddqM = ddqM.*dq;
    %
    %     ddqP = (12*dqMax.*tempExp2)./((tempExp2 + 1).^2.*(c1 - qMax));
    %     ddqP = ddqP.*dq;
    
    %% Conventional method
    %     dqM = max(kappa*(qMin - q), dqMin);
    %     dqP = min(kappa*(qMax - q), dqMax);
    
    jointLimit(index, :) = [dqP; dqM];
end

% figure;
% plot(t, jointLimit*180/pi);
% hold on;
% plot(t, jointVelocity(:, 4)*180/pi);
% 
% figure;
% plot(t, jointLimit(:, 11)*180/pi);
% hold on;
% plot(t, jointAngle(:, 4)*180/pi);
% 
% figure;
% hold on;
% plot(t, jointAngle(:, 4)*180/pi);
% plot(t, qMax(4)*ones(size(t))*180/pi);
% plot(t, qMin(4)*ones(size(t))*180/pi);
% 
% % for i = 1:length(t)
% %     theta = jointAngle(i, :);
% %     dtheta = jointVelocity(i, :);
% %     objective(i) = 0.5*dtheta*transpose(dtheta) + lambda*(theta - transpose(q0))*transpose(dtheta);
% % end
% % figure;
% % plot(t, objective);
% 
% % Distances: left, bottom, right, top
% figure;
% Rect = [0.075, 0.075, 0.9, 0.9];
% AxisPos = myPlotPos(1, 2, Rect);
% axes('Position', AxisPos(1, :));
% hold on;
% box on;
% h1 = line_fewer_markers(t, jointVelocity(:, 1), 6, '--r', 'MarkerFaceColor', 'g', 'MarkerSize', 6, 'LineWidth', 1);
% h2 = line_fewer_markers(t, jointVelocity(:, 2), 6, '-.b', 'MarkerFaceColor', 'g', 'MarkerSize', 6, 'LineWidth', 1);
% h3 = line_fewer_markers(t, jointVelocity(:, 3), 6, ':m', 'MarkerFaceColor', 'g', 'MarkerSize', 6, 'LineWidth', 2);
% ylim([-0.8, 0.8]);
% set(gca, 'FontSize', 14)
% legend(gca, {'data1','data2','data3'})
% a=axes('position',get(gca,'position'),'visible','off');
% hold on;
% box on;
% h4 = line_fewer_markers(t, jointVelocity(:, 4), 6, '-k', 'MarkerFaceColor', 'g', 'MarkerSize', 6, 'LineWidth', 1);
% h5 = line_fewer_markers(t, jointVelocity(:, 5), 6, '--rs', 'MarkerFaceColor', 'g', 'MarkerSize', 6, 'LineWidth', 1);
% h6 = line_fewer_markers(t, jointVelocity(:, 6), 6, '-.bs', 'MarkerFaceColor', 'g', 'MarkerSize', 6, 'LineWidth', 1);
% h7 = line_fewer_markers(t, jointVelocity(:, 7), 6, ':ms', 'MarkerFaceColor', 'g', 'MarkerSize', 6, 'LineWidth', 1);
% ylim([-0.8, 0.8]);
% set(gca, 'FontSize', 14);
% legend(a, {'data4','data5','data6','data7'})
% 
% axes('Position', AxisPos(2, :));
% hold on;
% box on;
% h8 = line_fewer_markers(t, jointVelocity(:, 4), 6, '-k', 'MarkerFaceColor', 'g', 'MarkerSize', 6, 'LineWidth', 1);
% h9 = line_fewer_markers(t, jointLimit(:, 4), 6, '-.bp', 'MarkerFaceColor', 'g', 'MarkerSize', 6, 'LineWidth', 1);
% h10 = line_fewer_markers(t, jointLimit(:, 11), 6, ':rp', 'MarkerFaceColor', 'g', 'MarkerSize', 6, 'LineWidth', 2);
% ylim([-1.5, 1.5]);
% set(gca, 'FontSize', 14);
% legend(gca, {'datax','datay','dataz'})
% set(gcf,'renderer','Painters')
% saveas(gcf,'velocityNCPZNN','epsc')

% save dataNCP.mat jointVelocity jointLimit t
Err_trace(:,5001)
jointVelocity(5001,:)
