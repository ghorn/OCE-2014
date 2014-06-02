function [f, grad_f, X] = fobj(U)

% get the dimension from the control vector
N = length(U);

% initial conditions
px0 = 0;
pz0 = -20;
vx0 = 10;
vz0 = 0;

% state history (returned to the user for plotting purposes)
X = zeros(4,N+1);
X(:,1) = [px0; pz0; vx0; vz0];

% simulate the system
for k=1:N
    [xnext,~,~] = integrate_airplane_ode(X(:,k), U(k));
    X(:,k+1) = xnext;
end

% objective function is final z velocity
f = X(4,end);

%% compute reverse-mode AD gradient
% We initialize the reverse loop, setting all bar quantities to zero, apart
% from the one corresponding to the selected output (already done for you) 

XBar = zeros(4,N+1);
UBar = zeros(1,N);
XBar(:,N+1) = [0;0;0;1];

% Now you need to write a reverse loop that computes XBar and UBar.
%
%
%
%
%

grad_f = UBar';

end
