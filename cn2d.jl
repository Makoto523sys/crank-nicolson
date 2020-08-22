using LinearAlgebra;
using Printf;

function linspace(start::Number, stop::Number, n)
	linx = zeros(Float64, n);
	delta = (stop - start) / (n - 1);
	for i = 1:n
		linx[i] = start + (i - 1) * delta;
	end
	return linx;
end

function tdma(n, a, b, c, r)
	ans = zeros(Float64, n);
	g = zeros(Float64, n);
	g[1] = b[1];
	s = zeros(Float64, n);
	s[1] = r[1];
	
	ans[1] = s[1] / g[1];
	
	for i = 2:n-1
		g[i] = b[i] - a[i] * c[i - 1] / g[i - 1];
		s[i] = r[i] - a[i] * s[i - 1] / g[i - 1];
		ans[i] = s[i] / g[i];
	end
	for i = n-2:1
		ans[i] = (s[i] - c[i] * ans[i + 1]) / g[i];
	end
	return ans
end

function init()
	global nx = 5;
	global ny = 5;
	global lx = Float64(100.0);
	global ly = Float64(100.0);
	global x = linspace(0, lx, nx);
	global y = linspace(0, ly, ny);
	global dx = lx / (nx - 1);
	global dy = ly / (ny - 1);
	global t = 0.0;
	global tlims = 300.0;;
	global dt = 10.0;
	global kappa = 0.835;
	global beta_x = kappa * dt / 2.0dx^2;
	global beta_y = kappa * dt / 2.0dy^2;
	global u = zeros(Float64, (ny, nx));
	for j in 1:ny
		u[j, 1] = 75.0;
		u[j, end] = 50.0;
	end
	for i in 1:nx
		u[1, i] = 0.0;
		u[end, i] = 100.0;
	end
	u[begin, begin] = (75.0 + 0.0) / 2.0;
	u[end, begin] = (75.0 + 100.0) / 2.0;
	u[begin, end] = (0.0 + 50.0) / 2.0;
	u[end, end] = (100.0 + 50.0) / 2.0;
end

function output(x::Array, y::Array, u::Array, t::Float64)
	filename = @sprintf("adi_cn2d_t=%03d.dat", t);
	open(filename, "w") do io
		for i in 1:ny
			for j in 1:nx
				write(io, @sprintf("%.5f %.5f %.5f\n", x[j], y[i], u[i, j]));
			end
		end
	end
end

function cn_direction_x(u::Array)
	a = zeros(Float64, nx - 2);
	b = zeros(Float64, nx - 2);
	c = zeros(Float64, nx - 2);
	r = zeros(Float64, nx - 2);
	for j = 1:ny-2
		for i = 1:nx-2
			a[i] = -beta_y;
			b[i] = 1.0 + 2.0beta_y;
			c[i] = -beta_y;
			r[i] = u[j+1, i+1] + beta_x*(u[j+2, i+1] - 2u[j+1, i+1] + u[j, i+1]);
		end
		r[begin] = r[begin] - a[begin] * u[j+1, begin];
		r[end] = r[end] - c[end] * u[j+1, end];
		u[j+1, begin+1:end-1] = tdma(nx - 2, a, b, c, r)
	end
end

function cn_direction_y(u::Array)
	a = zeros(Float64, nx - 2);
	b = zeros(Float64, nx - 2);
	c = zeros(Float64, nx - 2);
	r = zeros(Float64, nx - 2);
	for i in 1:nx-2
		for j in 1:ny-2
			a[j] = -beta_x;
			b[j] = 1.0 + 2.0beta_x;
			c[j] = -beta_x;
			r[j] = u[j+1, i+1] + beta_y * (u[j+1, i+2] - 2.0u[j+1, i+1] + u[j+1, i])
		end
		r[begin] = r[begin] - a[begin]*u[begin, i+1];
		r[end] = r[end] - c[begin]*u[end, i+1];
		u[begin+1:end-1, i] = tdma(ny-2, a, b, c, r);
	end
end

function main()
	init();
	@printf("nx = %d, ny = %d, dt = %.4f, dx = %.4f, dy = %.4f\n", nx, ny, dt, dx, dy);
	loop_count = 1;
	while t <= tlims
		@printf("t = %.5f\n", t);
		## 奇数回目はx方向を陰的にする.
		if loop_count % 2 == 1
			## x方向にcrank nicolson
			cn_direction_x(u);
			## y方向にcrank nicolson
			cn_direction_y(u);
			output(x, y, u, t);
		## 偶数回目はy方向を陰的にする.
		else
			## y方向にcrank nicolson
			cn_direction_y(u);
			## x方向にcrank nicolson
			cn_direction_x(u);
			output(x, y, u, t);
		end
		loop_count += 1;
		global t += dt;
	end
end

main()
