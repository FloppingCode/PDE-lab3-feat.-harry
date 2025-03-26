N = 100;
dx = 1/(N-2);
x = linspace(0,1,N-1);
f = @(t) 1;

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
        g = g(i);
    end
end



B = @(b) arrayfun(@(t) bj(b,t),x); 
dB = @(i) arrayfun(@(t) grad(t,i),x);
F = arrayfun(@(t) f(t), x);

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
delta = 0.05;
L = @(U,A,b) (transpose(F)*(U+A)-transpose(B(b)*D*U)*D*A);
% B is symmetric
%grad_L = transpose(DU)*dB(b)*D*A*dx;



b = [0.5 0.4];

err = 1;
tol = 1e-4;
while err > tol
    W = D\F'; 

    U = (diag(B(b))*D)\W;

    grad_L1 = (transpose(D*U)*diag(dB(1)))*D*U*dx;
    grad_L2 = (transpose(D*U)*diag(dB(2)))*D*U*dx;
    disp(["gradients", grad_L1,grad_L2, "b", b]);

    bnxt = b - delta*[grad_L1, grad_L2];
    err = norm(bnxt - b);
    b = bnxt;
    %disp(b);
end

disp(b);
