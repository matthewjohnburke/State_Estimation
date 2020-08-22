%% Extended Kalman Filter
%{ 
The following is an implementation of an Extended Kalman Filter. To demonstrate
its performance, the following code simulates the rotational movement of a satellite. The dynamics of the vehicle are
The satellite was modeled using Euler's equations of rigid body motion. 
The EKF provides estimates of the angular velocities of the nonlinear system.
This code was created as a course assignemt for AA273: State Estimation and Filtering for Aerospace
Applications.
%}

clc;
clear;

% Define Time Vectors
t_f = 15;
delta_t = 0.001;
time = 1:delta_t:t_f;
N = length(time);

% Define Parameters
Jx = 1;
Jy = 5;
Jz = 5;
J = [Jx, 0, 0; 0, Jy, 0; 0, 0, Jz];
Tau = [0;0;0];
c = 10;
Q = 0.004 * eye(3);
R = 0.1*eye(3);
mu_V = [0;0;0];
mu_W = [0;0;0];

% Initialize Space
wx = zeros(1,N);
wy = zeros(1,N);
wz = zeros(1,N);
Sigma = zeros(3,3,N);
mu = zeros(3,N);
A = zeros(3,3,N);
W = zeros(3,N);
V = zeros(3,N);
C = zeros(3,3,N);
Y = zeros(3,N);
Y_h = Y;

% Set Initial Conditionsd
wx(1) = 10;
wy(1) = 0.1;
wz(1) = 0.1;
mu(:,1) = [10,0,0]';
Sigma(:,:,1) = diag([1 1 1]);
W(:,1) = mvnrnd(mu_V,Q)';
V(:,1) = mvnrnd(mu_V,R)';

A(:,:,1) = [1, delta_t*(J(2,2) - J(3,3)) * mu(1,1), delta_t*(J(2,2) - J(3,3))*mu(2,1);...
         delta_t*(J(3,3) - J(1,1))*mu(3,1), 1, delta_t*(J(3,3) - J(1,1))*mu(1,1);...
         delta_t*(J(1,1) - J(2,2))*mu(2,1), delta_t*(J(1,1) - J(2,2))*mu(1,1), 1];


% Define States
X = [wx;wy;wz];

for j = 1:3
    if abs(X(j,1)) < c
        C(j,j,1) = 1;
        Y(j,1) = X(j,1);
        Y_h(j,1) = mu(j,1);
    elseif abs(X(j,1)) >= c
        C(j,j,1) = 0;
        Y(j,1) = sign(X(j,1)) * c;
        Y_h(j,1) = sign(mu(j,1))*c;
    end
end

Y(:,1) = Y(:,1) + V(:,1);
Y_h(:,1) = Y_h(:,1) + V(:,1);

for t = 2:(N)
   
    % EKF Predict
    mu(:,t) = [delta_t*(((J(2,2) - J(3,3))/J(1,1)) * (mu(2,t-1)*mu(3,t-1)) + Tau(1,1)) + mu(1,t-1);...
        delta_t*(((J(3,3) - J(1,1))/J(2,2)) * (mu(3,t-1)*mu(1,t-1)) + Tau(2,1)) + mu(2,t-1);...
        delta_t*(((J(1,1) - J(2,2))/J(3,3)) * (mu(1,t-1)*mu(2,t-1)) + Tau(3,1)) + mu(3,t-1)];
    Sigma(:,:,t) = A(:,:,t-1)*Sigma(:,:,t-1)*A(:,:,t-1)' + Q;
    
    % State and Measurement Calculations
    V(:,t) = mvnrnd(mu_V,R)';
    W(:,t) = mvnrnd(mu_V,Q)';
    
    X(:,t) = [delta_t*(((J(2,2) - J(3,3))/J(1,1)) * (X(2,t-1)*X(3,t-1)) + Tau(1,1)) + X(1,t-1);...
        delta_t*(((J(3,3) - J(1,1))/J(2,2)) * (X(3,t-1)*X(1,t-1)) + Tau(2,1)) + X(2,t-1);...
        delta_t*(((J(1,1) - J(2,2))/J(3,3)) * (X(1,t-1)*X(2,t-1)) + Tau(3,1)) + X(3,t-1)] + W(:,t);
    
    for j = 1:3
        if abs(X(j,t)) < c
            Y(j,t) = X(j,t);
        elseif abs(X(j,t)) >= c
            Y(j,t) = sign(X(j,t)) * c;
        end
    end
    
    for j = 1:3
        if abs(mu(j,t)) < c
            Y_h(j,t) = mu(j,t);
        elseif abs(mu(j,t)) >= c
            Y_h(j,t) = sign(mu(j,t))*c;
        end
    end
    Y(:,t) = Y(:,t) + V(:,t);
    
    

    for j = 1:3
        if abs(mu(j,t)) < c
            C(j,j,t) = 1;
        elseif abs(mu(j,t)) >= c
            C(j,j,t) = 0;
        end
    end    
    
    % EKF Update
    mu(:,t) = mu(:,t) + Sigma(:,:,t) * C(:,:,t)'*inv(C(:,:,t)*Sigma(:,:,t)*C(:,:,t)' + R) * (Y(:,t) - Y_h(:,t));
    Sigma(:,:,t) = Sigma(:,:,t) - Sigma(:,:,t)*C(:,:,t)'*inv(C(:,:,t)*Sigma(:,:,t-1)*C(:,:,t)' + R) * C(:,:,t) * Sigma(:,:,t);
    
    % State A Matrix Calculations
    A(:,:,t) = [1, delta_t*(J(2,2) - J(3,3)) * mu(3,t), delta_t*(J(2,2) - J(3,3))*mu(2,t);...
         delta_t*(J(3,3) - J(1,1))*mu(3,t), 1, delta_t*(J(3,3) - J(1,1))*mu(1,t);...
         delta_t*(J(1,1) - J(2,2))*mu(3,t), delta_t*(J(1,1) - J(2,2))*mu(1,t), 1];
    
end

figure(1);                                          
subplot(3,1,1); grid on;
    plot(time,X(1,:),time,mu(1,:))
    title('State vs Expected Values'); 
    xlabel('time'); ylabel('Angular Velocity');
    legend('wx','Expected wx')
subplot(3,1,2); grid on
    plot(time,X(2,:),time,mu(2,:))
    xlabel('time'); ylabel('Angular Velocity');
    legend('wy','Expected wy')
subplot(3,1,3); grid on
    plot(time,X(3,:),time,mu(3,:))
    xlabel('time'); ylabel('Angular Velocity');
    legend('wz','Expected wz')
