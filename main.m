N = 100;
dx = 1/(N);
x = linspace(dx, 1-dx, N-1);
%f = @(t) 1;
intf = @(x) x^2/2; 

function x = bj(b,t)
	if t < (1/3)
		x = b(1);
    elseif t <= (2/3) 
		x = b(2);	
	else
		x = (1-b(1)-b(2));
	end
end

function g = grad(t,i) 
	if t < (1/3)
		g = [1, 0];
        g = g(i);
    elseif t <= (2/3)
		g = [0, 1];
        g = g(i);
	else
		g = [-1,-1];
        %g = [0,0];
        g = g(i);
    end
end



B = @(b) arrayfun(@(t) bj(b,t),x); 
dB = @(i) arrayfun(@(t) grad(t,i),x);
F = arrayfun(@(t) intf(t), x);

%F(1) = 0;
%F(end) = 0;

% D matrix
D = zeros(N-1,N-1);
for i = 1:N-1
	D(i,i) = -2;
end	
for i = 1:N-2
    D(i,i+1) = 1;
	D(i+1,i) = 1;	
end
D = D*(1/dx)^2;	
% how we solve§
% W = D\F;
% U = (B(b)*D)\W;
% bestäm b_i nummeriskt
% bestäm B från (b1,b2)
delta = 0.5;
L = @(U,A,b) (transpose(F)*(U+A)-transpose(B(b)*D*U)*D*A);
% B is symmetric
%grad_L = transpose(DU)*dB(b)*D*A*dx;


it = 0;
iters = [];
total = [];

b = [0.2 0.2];

bb = B(b);

err = 1;
tol = 1e-6;
while err > tol
    W = D\F'; 


    U = (diag(B(b))*D)\W;

    iters(end+1) = it;
    total(end+1) = (F)*U*dx;
    it = it + 1;
    %U(1) = 0;
    %U(end) = 0;

    grad_L1 = -(transpose(D*U)*diag(dB(1)))*D*U*dx;
    grad_L2 = -(transpose(D*U)*diag(dB(2)))*D*U*dx;
    disp(["gradients", grad_L1,grad_L2, "b", b]);

    bnxt = b - delta*[grad_L1, grad_L2];

    %bnxt(1) = max(bnxt(1), 1e-3);  % b_1 > 0
    %bnxt(2) = max(bnxt(2), 1e-3);  % b_1 > 0
    %if sum(bnxt) >= 1
    %    bnxt = bnxt / (sum(bnxt) + 1e-3);  % scale to keep sum < 1
    %end

    err = norm(bnxt - b);
    b = bnxt;
    %disp(b);
end

disp(b);

disp(1-b(1)-b(2));

plot(iters, total);
xlabel("Iteration");
ylabel("F^T*U*dx");
%title("g = 1, a(x) = 1+x, f(x) = 0")
%legend("u(x) Temperature","exact solution");
