#include <cmath>
#include <iostream>
#include <vector>
//#include <bessel.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_math.h>
#include "matplotlibcpp.h"
#include <Python.h>
#include <opencv2/opencv.hpp>

namespace plt = matplotlibcpp;

#define X_MAX 20
#define dx 1e-3
#define g 9.8

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

double get_u(double w, double x, double l, double A, double t, std::vector<double>& mkv)
{
    double numeratorFirst = gsl_sf_bessel_J0(2 * w * std::sqrt(x / g));
    double denominatorFirst = gsl_sf_bessel_J0(2 * w * std::sqrt(l / g));
    double termFirst = A / (w * w) * std::sin(w * t) * (numeratorFirst / denominatorFirst - 1);
    double sum = 0;
    for (int i = 0; i < mkv.size(); i++)
    {
        double wk = mkv[i] * std::sqrt(g) / (2 * std::sqrt(l));
        double numeratorSecond = gsl_sf_bessel_J0(mkv[i] * std::sqrt(x / l)) * std::sin(wk * t);
        double denominatorSecond = (wk * wk - w * w) * mkv[i] * mkv[i] * gsl_sf_bessel_J1(mkv[i]);
        sum += numeratorSecond / denominatorSecond;
    }
    double termSecond = -4 * A * w * std::sqrt(l / g) * sum;
    double result = termFirst + termSecond;
    return result;
}

std::vector<double> get_expansion(double w, std::vector<double>& xv, double l, double A, double t, std::vector<double>& mkv)
{
    std::vector<double> uv;
    for (double i = 0; i < l / dx; i++)
    {
        uv.push_back(get_u(w, xv[i], l, A, t, mkv));
    }
    return uv;
}

std::pair<std::vector<double>, std::vector<double>> get_end(std::vector<double>& xvl, std::vector<double>& y)
{
    double y2, y1, dy, dl;
    double sum = 0;
    int index;
    int n = y.size();
    for ( int k = 1; k < n - 1; k++ )
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

std::vector<double> get_xvl(double l)
{
    std::vector<double> xvl;
    for(double i = 0; i < l; i+=dx)
    {
        xvl.push_back(i);
    }
    return xvl;
}

// forcedOsc.exe w1 T1
int main(int argc, char *argv[])
{
    double arg_w1 = atof(argv[1]);
    double arg_T1 = atof(argv[2]);
    //std::cout << arg_w1 << std::endl;
    //std::cout << arg_T1 << std::endl;
    std::vector<double> temp;
    std::cout << temp.max_size() << std::endl;

    double l = 1;
    double A = 1;
    std::vector<double> xv = get_xv();
    std::vector<double> yv = get_yv(xv);
    std::vector<double> mkv = get_mkv(xv.data(), yv.data());
    std::vector<double> xvl = get_xvl(l);

    double T1 = 4 * M_PI * std::sqrt(l) / (std::sqrt(g) * mkv[0]);   
    double w = arg_w1; // 1; // TODO: Arg
    double f = w / (2 * M_PI);
    double T = arg_T1 / f; // 1 / f; // TODO: Arg
    double k = 1;

    plt::figure();
    int j = 1;
    for (double t = 0; t <= T; t+=T/7)
    {
        std::vector<double> y = get_expansion(w, xv, l, A, t, mkv);       
        auto [xn, yn] = get_end(xvl, y);
        plt::clf();
        plt::plot(yn, xn);
        plt::ylim(-0.1, 1.1);
        plt::xlim(-0.2, 0.2);
        std::string title = "t = " + std::to_string(t);
        plt::title(title);
        plt::grid(true); 
        //plt::show();
        //std::cout << "First stage: " << t << " out of " << T << std::endl;
        std::string filename = "forcedFigures/plot_" + std::to_string(j) + ".png";
        plt::savefig(filename);

        j++;
    }   

    plt::figure();
    w = 2 * M_PI; // TODO: Arg
    T = 1;
    j = 1;
    for (double t = 0; t <= T; t+=T/7)
    {
        std::vector<double> y = get_expansion(w, xv, l, A, t, mkv);       
        auto [xn, yn] = get_end(xvl, y);
        plt::clf();
        plt::plot(yn, xn);
        plt::ylim(-0.1, 1.1);
        plt::xlim(-0.2, 0.2);
        std::string title = "t = " + std::to_string(t);
        plt::title(title);
        plt::grid(true); 
        //plt::show();
        //std::cout << "Second stage: " << t << " out of " << T << std::endl;
        std::string filename = "forcedFigures/plot_w2PI_" + std::to_string(j) + ".png";
        plt::savefig(filename);
        j++;
    }   

    plt::figure();
    double w1 = std::sqrt(g) * mkv[0] / (std::sqrt(l) * 2);
    w = w1 - w1 * 0.1;
    T = 1; //10;
    j = 1;
    for (double t = 0; t <= T; t+=T/7)
    {
        std::vector<double> y = get_expansion(w, xv, l, A, t, mkv);       
        auto [xn, yn] = get_end(xvl, y);
        plt::clf();
        plt::plot(yn, xn);
        plt::ylim(-0.1, 1.1);
        plt::xlim(-1., 1.);
        std::string title = "t = " + std::to_string(t);
        plt::title(title);
        plt::grid(true); 
        //plt::show();
        //std::cout << "Third stage: " << t << " out of " << T << std::endl;
        std::string filename = "forcedFigures/plot_closeToResonance_" + std::to_string(j) + ".png";
        plt::savefig(filename);
        j++;
    }  
}