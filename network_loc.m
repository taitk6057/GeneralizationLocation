
function [pos_free] = network_loc(N, E, pos_anchor, rho);
N=50; R=0.4; s=0.05;
[E, pos, K] = network_loc_data(N, R);
pos_anchor = pos(N-K+1:N, :);
L = size(E,1);
d = sqrt(sum( (pos(E(:,1),:) - pos(E(:,2),:)).^2, 2));
rho = (1 + s*randn(L,1)) .* d;

lambda = 1.0;
u = rand(N-K,1);
v = rand(N-K,1);
theta = [u;v];
cost=[];

for iters = 1:100
    
    f = (theta.^2);
	Df = (2*theta);
    cost = [cost, norm(2*Df'.*f).^2];
	
	if norm(2*Df'.*f) <= 10e-5; break; end;
    
	newtheta = theta - [Df'; sqrt(lambda)*eye(82)] \ [f;zeros(1,1)];
    newf = newtheta.^2;

	if norm(newf) < norm(f)
		theta = newtheta;
		lambda = .8 * lambda; 
	else
		lambda = 2 * lambda;
	end;
end;

%plot
semilogy([0:iters-1], cost, 'o', [0:iters-1], cost, '--');
xlabel('k');
ylabel('cost function');
grid on;