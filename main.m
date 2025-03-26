N = 100
dx = 1/N
x = linespace(0,1,dx)
f = @(t) 1

function bj(b,t)
	if t < (1/3)
		return b(1);
	elseif t <= (2/3) 
		return b(2);	
	else
		return (1-b(1)-b(2));
	end
end

function grad(b,t) 
	if t < (1/3)
		return 1, 0;
	elseif t <= (2/3)
		return 0, 1;
	else
		return -1,-1;
end
B = @(b) arrayfun(@(t) bj(b,t),x); 
dB = @(b) arrayfun(@(t) grad(b,t),x);
F = arrayfun(@(t) f(t), x);
% D matrix
D = zeros(N-1,N-1);
for i = 1:N-1
	D(i,i) = -2;
	D(i,i+1) = 1;
	D(i+1,i) = 1;	
end	
D = D*(1/dx)^2	
% how we solve§
% W = D\F;
% U = (B(b)*D)\W;
% bestäm b_i nummeriskt
% bestäm B från (b1,b2)
delta = 0.5;
L = @(U,A,b) (transpose(F)*(U+A)-transpose(B(b)*D*U)*D*A);
% B is symmetric
grad_L = transpose(DU)*dB(b)*D*A*dx;
