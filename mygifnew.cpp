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
#define dx 1e-3

std::vector<double> get_xv()
{
    std::vector<double> xv;
    for(double i = 0; i < X_MAX; i+=dx)
    {
        xv.push_back(i);
    }
    return xv;
}

std::vector<double> get_yv(std::vector<double> xv) // how to do *xv
{
    std::vector<double> yv;
    for(double i = 0; i < X_MAX / dx; i++)
    {
        yv.push_back(gsl_sf_bessel_J0(xv[i]));
    }
    return yv;
}

std::vector<double> get_mkv(double* xv, double *yv)
{
    std::vector<double> mkv;
    //xv_count = xv.size();
    for (int k = 0; k < X_MAX / dx - 1; k++)
    {
        if ( yv[k] * yv[k + 1] <= 0)
        {
            mkv.push_back(xv[k]);
        }
    }
    return mkv;
}

double get_f(double k_coef, double b, double x)
{
    return k_coef * x + b;
}

// std::vector<double> get_akv(double* xv, double l, std::vector<double> mkv, double k_coef, double b )
// {
//     std::vector<double> akv;
//     double result = 0;
//     for (int k = 0; k < mkv.size(); k++ )
//     {
//         result = 0;
//         for (int i = 0; i < l / dx; i++)
//         {
//             result += get_f(k_coef, b, xv[i]) * gsl_sf_bessel_J0(mkv[k] * std::sqrt(xv[i] / l)) * dx;
//         }
//         result /= std::pow(l * gsl_sf_bessel_J1(mkv[k]), 2);
//         akv.push_back(result);
//     }
//     return akv;
// }

std::vector<double> get_akv(double* xvl, double l, std::vector<double> mkv, double k_coef, double b )
{
    std::vector<double> akv;
    double result = 0;
    //std::cout << "akv:";
    for (int k = 0; k < mkv.size(); k++ )
    {
        result = 0;
        for (int i = 0; i < l / dx - dx; i++)
        {
            result += get_f(k_coef, b, xvl[i]) * gsl_sf_bessel_J0(mkv[k] * std::sqrt(xvl[i] / l)) * dx;
        }
        result /= (l * std::pow(gsl_sf_bessel_J1(mkv[k]), 2));
        akv.push_back(result);
        //std::cout << akv[k] << " ";
    }
    return akv;
}

double get_item(double x, double* akv, std::vector<double>& mkv, double l, double a, double t)
{
    double result = 0;
    for(int k = 0; k < mkv.size(); k++)
    {
        result += akv[k] * gsl_sf_bessel_J0(mkv[k] * std::sqrt(x / l)) * std::sin(a * mkv[k] * t / (2 * std::sqrt(l) + M_PI / 2));
    }
    return result;
}

std::vector<double>get_expansion(double* akv, double l, double *xv, std::vector<double>& mkv, double k_coef, double b, double t, double a)
{
    //std::vector<double> akv = get_akv(xv, l, mkv, k_coef, b);
    std::vector<double> uv;
    // std::cout << "k_coef = " << k_coef << std::endl;
    for(int i = 0; i < l / dx; i++)
    {
        uv.push_back(get_item(xv[i], akv, mkv, l, a, t));
    }
    return uv;
}

std::vector<double> get_xvl(double l)
{
    std::vector<double> xvl;
    for(double i = 0; i < l; i+=dx)
    {
        xvl.push_back(i);
    }
    return xvl;
}

std::pair<std::vector<double>, std::vector<double>> get_end(std::vector<double>& xvl, std::vector<double>& y)
{
    double y2, y1, dy, dl;
    double sum = 0;
    int index;
    int n = y.size();
    for ( int k = 1; k < n; k++ )
    {
        y2 = y[n - k + 1];
        y1 = y[n - k];
        dy = std::abs(y2 - y1);
        dl = std::sqrt((dy * dy) + (dx * dx));
        sum += dl;
        if (sum >= 1.0)
        {
            index = n - k;
            break;
        }
    }
    std::vector<double> xn(xvl.begin() + index, xvl.end());
    std::vector<double> yn(y.begin() + index, y.end());
    //std::cout << "Size of xn" << xn.size() << std::endl << "Size of yn" << yn.size() << std::endl;
    return std::make_pair(xn, yn);
}

int main()
{
    double k_coef = -0.36;
    double b = -k_coef;
    double x = 0.1;
    double l = 1.;
    double a = std::sqrt(9.8);
    //double t;
    double t = 0.5;
    std::vector<double> xv = get_xv();
    std::vector<double> yv = get_yv(xv);
    std::vector<double> mkv = get_mkv(xv.data(), yv.data());
    std::vector<double> akv = get_akv(xv.data(), l, mkv, k_coef, b);
    //double result = get_item(x, akv.data(), mkv, l, a, t);
    //std::vector<double> uv = get_expansion(akv.data(), l, xv.data(), mkv, k_coef, b, t, a);
    std::vector<double> xvl = get_xvl(l);
    //std::vector<double> uv = get_expansion(akv.data(), l, xvl.data(), mkv, k_coef, b, t, a);
    std::vector<double> f0;
    for (int i = 0; i < xvl.size(); i++)
    {
        f0.push_back(get_f(k_coef, b, xvl[i]));
    }
    double T = 4 * M_PI * std::sqrt(l) / (a * mkv[0]);
    
    // double initial_t = 0;
    // std::vector<double> y_initial = get_expansion(akv.data(), l, xv.data(), mkv, k_coef, b, initial_t, a);
    // auto [xn_initial, yn_initial] = get_end(xvl, y_initial);

    // Плот для начального состояния
    // plt::plot(yn_initial, xn_initial, "b"); // Синий график
    // plt::plot(f0, xvl, "r"); // Красный график
    // plt::ylim(-0.1, 1.1);
    // plt::xlim(-0.4, 0.4);
    // plt::title("Initial state (t = 0)");
    // plt::grid(true);
    // plt::show();
    for (double t = 0; t < T; t+=T/7)
    {
        std::vector<double> y = get_expansion(akv.data(), l, xv.data(), mkv, k_coef, b, t, a);
        auto [xn, yn] = get_end(xvl, y);
        // std::cout << "xn:" << std::endl;
        // for (int i = 0; i < 100; i++)
        // {
        //     std::cout << xn[i] << std::endl;
        // }
        // std::cout << "yn:" << std::endl;
        // for (int i = 0; i < 100; i++)
        // {
        //     std::cout << yn[i] << std::endl;
        // }
        //plt::plot(xv, yv);
        plt::plot(yn, xn);
        plt::plot(f0, xvl);
        plt::ylim(-0.1, 1.1);
        plt::xlim(-0.4, 0.4);
        std::string title = "t = " + std::to_string(t);
        plt::title(title);
        plt::grid(true); 
        plt::show();
    }
    //std::cout << "f = " << f << std::endl;
    //std::cout << "Size of xv: " << xv.size() << std::endl;
    //std::cout << "Size of yv: " << yv.size() << std::endl;
    //std::cout << "Size of mkv: " << mkv.size() << std::endl;
    //std::cout << "Size of akv: " << akv.size() << std::endl;
    //std::cout << "Result of get_item function:" << result;
    //std::cout << "Size of uv:" << uv.size();
    //std::cout << "Size of xvl:" << xvl.size() << std::endl;
    //std::cout << "Size of f0:" << f0.size() << std::endl;
    //std::cout << "T = " << T << std::endl;
    //std::cout << "Size of xn" << xn.size() << std::endl << "Size of yn" << yn.size() << std::endl;
}