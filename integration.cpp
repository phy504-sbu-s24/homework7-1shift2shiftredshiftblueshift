#include <cmath>
#include <functional>
#include <iostream>
#include <random>
#include <iomanip>

double EPSILON = std::pow(10,-12); //epsilon for the sin^2 integrand

double trapezoid(std::function<double(double)> f, double x_down, double x_up, int N);
double monte_carlo(std::function<double(double)> f, double x_down, double x_up, int N);
double gaussian(double x);
double sin2(double x);

double trapezoid(std::function<double(double)> f, double x_down, double x_up, int N) { //Function to integrate a function using trapezoid rule
    //f is the integrand, x_down is the lower limit, x_up is the upper limit, and N is the number of slabs to divide the domain into
    double DX = std::abs((x_up - x_down))/N; //slab width
    double sum = 0.0;
    for (double x = x_down; x <= x_up; x+=DX ) {
        sum += f(x) + f(x+DX);
    }
    return 0.5*DX*sum;
}

double monte_carlo(std::function<double(double)> f, double x_down, double x_up, int N) { //Function to integrate a function using monte-carlo
    //f is the integrand, x_down is the lower limit, x_up is the upper limit, and N is the number of samples to draw
    std::mt19937 generator(501); //Choose a seed for reproducible results
    double sum{0.0};
    double x; //
    std::uniform_real_distribution<double> uniform(x_down, x_up); //Distribution from which to draw locations at which to evaluate f
    for (int i=0; i < N; i++) {
        x = uniform(generator);
        sum += f(x);
    }

    return sum*(x_up - x_down)/N;

}

double gaussian(double x) {
    return std::exp(-x*x);
}

double sin2(double x) {
    return std::pow(std::sin(1/(x*(2-x)+EPSILON)),2);
}

void main() {
    int N_max = 4096;
    std::cout << "Gaussian integration table" <<std::endl;
    std::cout << std::setw(15) << "Num samples" << std::setw(15) <<  "Trapezoid" << std::setw(15) << "Monte Carlo" << std::endl;
    int N=8;
    while (N <= N_max) {
        std::cout <<std::setw(15) << N <<std::setw(15) << trapezoid(gaussian, -5, 5, N) <<std::setw(15) << monte_carlo(gaussian, -5, 5, N) <<std::endl;
        N *= 2;
    }

    std::cout << "Sin^2 function integration table" <<std::endl;
    std::cout << std::setw(15) << "Num samples" << std::setw(15) <<  "Trapezoid" << std::setw(15) << "Monte Carlo" << std::endl;
    N=8;
    while (N <= N_max) {
        std::cout <<std::setw(15) << N <<std::setw(15) << trapezoid(sin2, 0, 2, N) <<std::setw(15) << monte_carlo(sin2, 0, 2, N) <<std::endl;
        N *= 2;
    }
}