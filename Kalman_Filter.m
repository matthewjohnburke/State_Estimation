%% Kalman Filter
%{ 
The following is an implementation of a basic Kalman Filter. To demonstrate
its performance, the following code simulates the movement of a 2D vehicle
where velocity and position are controlled. The dynamics of the vehicle are
simulated with noise defined below and the Kalman Filter provides estimates
for the position and velocity of the vehicle. This code was created as a
course assignemt for AA273: State Estimation and Filtering for Aerospace
Applications.
%}

delta_t = 1; % Time step in seconds
t_f = 20;    % Final time
time = 1:delta_t:t_f;

Q = eye(2);  % Covariance of W noise
mu_W = [0;0];
R = 9*eye(2); % Covariance of V noise
mu_V = [0;0];

P_current = [1000; 0]; % Initial Position in meters
S_current = [0; 50];   % Initial Velocity in m/s

Sigma(:,:,1) = [250000*eye(2),zeros(2); zeros(2), eye(2)];
mu(:,1) = [1500;100;0;55];


X(:,1) = [P_current;S_current];

A = [1,0,delta_t,0;0,1,0,delta_t;0,0,1,0;0,0,0,1];
B = [0,0;0,0;1,0;0,1];
C = [1, 0, 0, 0; 0,1,0,0];

W(:,1) = mvnrnd(mu_W,Q);
V(:,1) = mvnrnd(mu_V,R);
Y(:,1) = C * X(:,1) + V(:,1);

for t = 1:delta_t:(t_f-1)
    u(:,t) = -2.5*[cos(0.05*t); sin(0.05*t)];
    
    
    % KF Update
    mu(:,t) = mu(:,t) + Sigma(:,:,t) * C' * inv(C * Sigma(:,:,t) * C' + R) * (Y(:,t) - C*mu(:,t));
    Sigma(:,:,t) = Sigma(:,:,t) - Sigma(:,:,t) * C' * inv(C*Sigma(:,:,t)*C' + R)*C*Sigma(:,:,t);
   
    
    figure(1);

    [ellipseconst,w] = Ellipse_plot( mu(1:2,t), Sigma(1:2,1:2,t) );
    [ellipseconst,w] = Ellipse_plot_v( X(1:2,t),mu(3:4,t), Sigma(3:4,3:4,t) );
    
    W(:,t) = mvnrnd(mu_W,Q);
    V(:,t) = mvnrnd(mu_V,R);

    X(:,t+1) = A * X(:,t) + B * u(:,t) + [0;0; W(:,t)]; % State Update
    Y(:,t+1) = C * X(:,t) + V(:,t);                    % Helicopter position measurements
    
    % KF Predict
    mu(:,t+1) = A * mu(:,t) + B * u(:,t);
    Sigma(:,:,t+1) = A * Sigma(:,:,t) * A' + [zeros(2),zeros(2);zeros(2),Q];
    
end

figure(1)
plot(X(1,:),X(2,:),'g')
xlabel('X [m]'); ylabel('Y [m]');
legend('Position Error Ellipses','Velocity Error Ellipses','Helicopter Position');

%% Ellipse Plotting Functions


function [ellipseconst,w] = Ellipse_plot( mu, Sigma )

P = 0.95;
ellipseconst = -2 * log((1-P)*sqrt(det(Sigma)));
r  =  sqrt(abs(ellipseconst));
theta =0:0.1:2*pi;
w1 = r * cos(theta);
w2 = r * sin(theta);
w = [w1; w2];
x = sqrtm(Sigma) * w  +  mu*ones(1,size(w,2));
plot(x(1,:),x(2,:),'g'); hold on;

end

function [ellipseconst,w] = Ellipse_plot_v( x, mu, Sigma )

P = 0.95;
ellipseconst = -2 * log((1-P)*sqrt(det(Sigma)));
r  =  sqrt(abs(ellipseconst));
theta =0:0.1:2*pi;
w1 = r * cos(theta);
w2 = r * sin(theta);
w = [w1; w2];
x = sqrtm(Sigma) * w  +  mu*ones(1,size(w,2)) + x;
plot(x(1,:),x(2,:),'b'); hold on;

end

