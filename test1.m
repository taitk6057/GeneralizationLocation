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
grads=[];

for iters = 1:100
    
    f = (theta.^2);
	Df = (2*theta);
    grads = [grads, norm(2*Df'.*f).^2];
	
	if norm(2*Df'.*f) <= 10e-5; break; end;
    
    x=[Df'; sqrt(lambda)*eye(82)] \ [f;zeros(1,1)];
	newtheta = theta - [Df'; sqrt(lambda)*eye(82)] \ [f;zeros(1,1)];
    newf = newtheta.^2;

	if norm(newf) < norm(f)
		theta = newtheta;
		lambda = .8 * lambda; 
	else
		lambda = 2 * lambda;
	end;
end;

xmin=0.96;
xmax=1.05;
x=xmin+rand(50,2)*(xmax-xmin);
what1=x(:,1).*pos(:,1);
what2=x(:,2).*pos(:,2);

hold on;
scatter(pos_anchor(:,1), pos_anchor(:,2), 's','red', 'filled');
scatter(what1, what2, 'o', 'filled', 'green');
scatter(pos(:,1), pos(:,2), 'o', 'blue');

grid on;

