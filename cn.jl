using LinearAlgebra;
using Printf;

function linspace(start::Number, stop::Number, n::Int64)
	linx = zeros(Float64, n);
	delta = (stop - start) / (n - 1);
	for i = 1:n
		linx[i] = start + delta * (i - 1)
	end
	return linx;
end

lx = 1.0;
t = 0.0
nx = 11;
dx = lx / (nx - 1);
dt = 1e-2;
tlim = 0.1;
kappa = 1.0
beta = dt * kappa / 2.0dx^2;
t_save = linspace(0.0, tlim, nx);
save_count = 1;

function solve_linear(n::Int64, a::Float64, b::Float64, c::Float64, r::Array)
	A = zeros(n, n);
	for i = 2:n+1
		for j = 2:n+1
			if(i-1 == j-1-1)
				A[i-1, j-1] = c;
			elseif (i-1==j-1)
				A[i-1,j-1] = b;
			elseif (i-1-1 == j-1)
				A[i-1,j-1] = a;
			end
		end
	end
	return A \ r;
end

function crank_nicolson(u::Array)
	a = Float64(0.0);
	b = Float64(0.0);
	c = Float64(0.0);
	r = zeros(Float64, nx - 2);
	for i = 2:nx-1
		a = -beta;
		b = 1.0 + 2.0beta;
		c = -beta;
		r[i-1] = u[i] + beta*(u[i+1] - 2.0u[i] + u[i-1]);
	end
	r[begin] = r[begin] - a*u[begin]
	r[end] = r[end] - c * u[end];
	u[2:end-1] = solve_linear(nx - 2, a, b, c, r) 
end

function output(x::Array, u::Array, t::Float64)
	filename = @sprintf("cn_t=%4.3f.dat", t);
	open(filename, "w") do io
		for i = 1:nx
			write(io, @sprintf("%.5f\t%.5f\n", x[i], u[i]));
		end
	end;
end

function main()
	@printf("beta = %.5f, dt = %.5f, dx = %.5f\n", beta, dt, dx);
	u = ones(Float64, nx);
	u[begin] = u[end] = 0.0;
	x = linspace(0, lx, nx);
	while t < tlim
		output(x, u, t);
		crank_nicolson(u);
		global t = t + dt;
	end
end

main()
