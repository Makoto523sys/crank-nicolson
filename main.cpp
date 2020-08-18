#include<iostream>
#include<fstream>
#include<vector>
#include<algorithm>
#include<string>

double lx = 1.0;
int nx = 11;
double dt = 1./100;
double dx = lx / (nx - 1);
double kappa = 1.0;
double beta = (dt * kappa) / (2.0 * dx * dx);
double t = 0.0;
double tlim = 0.1;

std::vector<double> linspace(double begin, double end, int nx){
	std::vector<double> lin;
	double delta = (end - begin) / (nx - 1);
	for(int i = 0; i < nx; i++) lin.push_back(begin + i * delta);
	return lin;
}

void output(std::vector<double> x, std::vector<double> u, double t){
	std::string filename = "crank-nicolson" + std::to_string(t) + ".dat";
	std::ofstream out{filename};

	for(int i = 0; i < nx; i++)
		out << x[i] << "\t\t" << u[i] << std::endl;
	out.close();
}

std::vector<double> tdma(int n, std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> r){
	std::vector<double> ans(n);
	std::vector<double> g(n); g[0] = b[0];
	std::vector<double> s(n); s[0] = r[0];
	ans[0] = s[0] / g[0];
	for(int i = 1; i < n; i++){
		g[i] = b[i] - a[i] * c[i - 1] / g[i - 1];
		s[i] = r[i] - a[i] * s[i - 1] / g[i - 1];
		ans[i] = s[i] / g[i];
	}
	for(int i = n - 2; i >= 0; i--) ans[i] = (s[i] - c[i] * ans[i + 1]) / g[i];
	return ans;
}

void crank_nicolson(std::vector<double>& u){
	std::vector<double> a(nx - 2);
	std::vector<double> b(nx - 2);
	std::vector<double> c(nx - 2);
	std::vector<double> r(nx - 2);
	std::vector<double> up1(nx - 2);
	std::vector<double> u_tmp(nx); std::copy(u.begin(), u.end(), u_tmp.begin());
	for(int i = 1; i < nx - 1;i ++){
		a[i - 1] = -beta;
		b[i - 1] = 1. + 2. * beta;
		c[i - 1] = -beta;
		r[i - 1] = u[i] + beta * (u[i+1] - 2.*u[i] + u[i-1]);
	}
	r[0] = r[0] -a[0] * u[0];
	r[nx-3] = r[nx - 3] - c[nx-3]*u[nx-1];
	up1 = tdma(nx-2, a, b, c, r);
	std::copy(up1.begin(), up1.end(), u_tmp.begin() + 1);
	std::copy(u_tmp.begin() + 1, u_tmp.end() - 1, u.begin() + 1);
}

int main(void){
	auto x = linspace(0.0, lx, nx);
	std::vector<double> u(nx); std::fill(u.begin(), u.end(), 1.0);
	u[0] = u[nx-1] = 0.0;
	while(t < tlim){
		t += dt;
		crank_nicolson(u);
		output(x, u, t);
	}
	return 0;
}
