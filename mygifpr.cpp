#include <cmath>
#include <iostream>
#include <vector>
//#include <bessel.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_math.h>
#include "matplotlibcpp.h"
#include <Python.h>

namespace plt = matplotlibcpp;

#define X_MAX 20
//#define dx 1e-3
//#define xv 0:dx:X_MAX

double fun_f(const double x, const double k, const double b)
{
    return k * x + b;
}

std::vector<double> get_ak(double* mkv, int mk_count, double dx, int xv_count, double *xv, double k, double b, double l, double x)
{
    double ak;
    std::vector<double> akv;
    for (int k = 0; k < mk_count; k++)
    {
        ak = 0;
        for (int i = 0; i < l; i++)
        {
            ak += fun_f(x, k, b) * gsl_sf_bessel_J0(mkv[k] * std::sqrt(xv[i]/ l)) * dx;
        }
        ak /= (l * std::pow(gsl_sf_bessel_J1(mkv[k]), 2));
        akv.push_back(ak);
    }
    return akv;
}

//double get_u(double x, double t, double *mkv, int mk_count, double a, double l, double dx, int xv_count, double b, double* xv, double k, double w)
double get_u(double x, double t, double *mkv, int mk_count, double a, double l, double dx, int xv_count, double b, double* xv, double k)
{
    double sum = 0;
    double lk;
    std::vector<double> akv;
    akv = get_ak(mkv, mk_count, dx, xv_count, xv, k, b, l, x);
    for (int i = 0; i < mk_count; i++) {
        //lk = sqrt(pow(mkv[k], 2) / (4*l) - pow((w/a), 2));
        sum += akv[i] * cos(a * mkv[i] * t / (2 * std::sqrt(l))) * gsl_sf_bessel_J0(mkv[i] * std::sqrt(x / l));
        //sum += akv[k] * cos(a * lk * t) * gsl_sf_bessel_J0(mkv[k] * sqrt(x / l));
    }
    return sum;
}

//std::vector<double> get_uv(double *xv, double t, double *mkv, int mk_count, double a, double l, double dx, int xv_count, double b, double k, double w)
std::vector<double> get_uv(double *xv, double t, double *mkv, int mk_count, double a, double l, double dx, int xv_count, double b, double k)
{
    double sum = 0;
    std::vector<double> uv;

    //for(int i = 0; i < xv_count; i++)
    double limit = l / dx + 1;
    //std::cout << "l / dx = " << limit << std::endl;
    for(int i = 0; i < limit; i++)
    {
        //std::cout << "in" << std::endl;
        //uv.push_back(get_u(xv[i], t, mkv, mk_count, a, l, dx, xv_count, b, xv, k, w));
        uv.push_back(get_u(xv[i], t, mkv, mk_count, a, l, dx, xv_count, b, xv, k));
        //std::cout << uv[i];
    }
    //std::cout << uv.size() << std::endl;
    return uv;
}

std::pair<std::vector<double>, std::vector<double>> get_end(const std::vector<double>& xv, const std::vector<double>& yv, double dx)
{
    double sum = 0;
    int index = 0;
    int n = yv.size();
    //std::cout << "n = " << n << std::endl;
    //for(int i = 0; i < yv.size() - 1; i++ )
    //for (int k = 1; k < n; k++)
    for(int k = n - 1; k > 0; k--)
    {
        double y1, y2, dy, dl;
        std::vector<double> xn, yn;
        // y1 = yv[i];
        // y2 = yv[i + 1];
        y2 = yv[n - k];
        y1 = yv[n - k - 1];
        dy = std::abs(y2 - y1);
        dl = std::sqrt(dy * dy + dx * dx);
        sum += dl;
        if (sum >= 1)
        {
            // index = i;
            //index = n - k - 1;
            index = k;
            break;
        }
    }
    //std::cout << "l = " << sum << std::endl;
    //std::vector<double> xn(xv.begin(), xv.begin() + index);
    //std::vector<double> yn(yv.begin(), yv.begin() + index);
    std::vector<double> xn(xv.begin() + index, xv.end());
    std::vector<double> yn(yv.begin() + index, yv.end());
    //std::vector<double> xn(xv.end() - index, xv.end());
    //std::vector<double> yn(yv.end() - index, yv.end());
    return std::make_pair(xn, yn);
}

int main()
{
    double l, k_coef, b, t, a, dx, T;
    l  = 1;
    dx = 0.001;
    //k_coef = -0.01;
    k_coef = -0.36;
    b = -k_coef;
    a = std::sqrt(9.8);
    int mk_count, xv_count;

    std::vector<double> xv, yv, mkv, uv, y, f0;

    for (double x = 0; x < X_MAX; x += dx)
    {
        xv.push_back(x);
    } 
    xv_count = xv.size();

    for (size_t i = 0; i < xv_count; i++)
    {
        double x = xv[i];
        double y = gsl_sf_bessel_J0(x);
        yv.push_back(y);
    }
    double yv_count = yv.size();

    for (size_t k = 0; k < xv_count; k++)
    {
        if (yv[k] * yv[k + 1] <= 0)
        {
            double xk = xv[k];
            mkv.push_back(xk);
        }
    }
    mk_count = mkv.size();

    //std::cout << mk_count << std::endl;
    int stop;
    stop = l / dx + 1;

    for(int i = 0; i < stop; i++)
    {
        f0.push_back(fun_f(xv[i], k_coef, b));
    }

    T = 4 * M_PI * std::sqrt(l) / (a * mkv[0]);
    
    // std::vector<double> xvn;
    // double limit = l / dx + 1;
    // for(int i = 0; i < limit; i+=dx)
    // {
    //     xvn.push_back(xv[i]);
    // }
    // std::cout << "xvn = " ; //<< xvn.size();

    std::vector<double> xvn;
    xvn.assign(xv.begin(), xv.begin() + stop);
    for (double t = 0; t < T + T / 7; t += T / 7)
    {
        y = get_uv(xvn.data(), t, mkv.data(), mk_count, a, l, dx, xv_count, b, k_coef);
        std::cout << "xvn: " << std::endl;
        for(int j = 0; j < 100; j++)
        {
            std::cout << xvn[j] << " ";
        }
        std::cout << std::endl << "y: " << std::endl;
        for(int j = 0; j < 100; j++)
        {
            std::cout << y[j] << " ";
        }
        std::cout << std::endl;

        auto [xn, yn] = get_end(xvn, y, dx);
        
        std::cout << "xn: " << std::endl;
        for(int j = 0; j < 100; j++)
        {
            std::cout << xn[j] << " ";
        }
        std::cout << std::endl << "yn: " << std::endl;
        
        for(int j = 0; j < 100; j++)
        {
            std::cout << yn[j] << " ";
        }
        //std::cout << "y = " << y.size() << std::endl << "xv = " << xv.size();
        //std::cout << std::endl << "xvn = " << xvn.size();
        //plt::plot(y, xvn);
        plt::plot(yn, xn);
        plt::plot(f0, xvn);
        plt::ylim(-0.1, 1.1);
        plt::xlim(-0.4, 0.4);
        std::string title = "t = " + std::to_string(t);
        plt::title(title);
        plt::grid(true); 
        plt::show();
    }
    return 0;
}   

// int main() {
//     plt::plot({1,3,2,4});
//     plt::show();
//     std::cout << "Hooray!";
// }