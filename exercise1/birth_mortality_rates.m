function [beta,mu] = birth_mortality_rates()

n = 100; % max age
ages = (1:n)';

% Gompertz–Makeham law of mortality
alpha_ = 0.0005202;
beta_ = 0.000007786;
gamma_ = 0.1116;
x = ages(1:n-1);
mu = (beta_*exp(gamma_*x) + alpha_) .* exp(-alpha_*x - beta_/gamma_*(exp(gamma_*x) - 1));
%mu = exp(-alpha_*x - beta_/gamma_*(exp(gamma_*x) - 1));
%mu(end) = 1;
x = ages(1:n);
sigma_ = 6;
beta = 2/(sigma_*sqrt(2*pi)).*exp(-0.5*(x-32).^2./(sigma_^2));
