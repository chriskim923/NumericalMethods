// Newton's Method
//
// This program uses the Newton's Method algorithm to determine the root
//
//*****************
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <stdio.h>
using namespace std;
//*****************

// Declare functions
double f1(double x);
double f2(double x);
double df1(double x);
double df2(double x);

// Declaration of output stream

ofstream data_out("output.txt"); // output data

// Assign functions
// (i) All roots are contained within [0.5 : 10] based on plot
double f1(double x) {
    return sin(x) + log(x) - 1;
}

double df1(double x) {
    return 1/x + cos(x);
}

// (ii) All roots are contained within [-5 : 5] based on plot
double f2(double x) {
    return exp(x) + pow(x,2) - 3*x - 2;
}

double df2(double x) {
    return exp(x) + 2*x - 3;
}

int main() {
    // Variable declarations
    double x, xi;               // x & current guess x
    double f, fi, df;           // f & current f(x) & df values
    int N = 20;                 // max number of iterations
    int i;                      // loop counter
    double f_left, f_right;     // function values for bracketing root
    double x_left, x_right;     // x values for the brackets
    double step;                // loop counter for bracketing root
    
    const float tolerance = 1e-6;
    
    // Root finder
    // Loop through initial range to find bracketed range of each root
    // Part (i)
    data_out<<"Roots found for part (i) "<<endl;
    for(step = 0.5; step < 10; step+=0.5){
        x_left = step;
        x_right = step+0.5;
        f_left = f1(x_left);
        f_right = f1(x_right);
        if(f_left * f_right < 0) {          // opposite signs
            xi = (x_left + x_right)/2;		// guess for Newton root finding: average between left and right brackets
            // Loop through Newton-Raphson
            for(i = 0; i < N; i++){
                x = xi;
                f = f1(x);
                df = df1(x);
                // Newton step
                xi = x - f/df;
                fi = f1(xi);
                // Check if convergence is within tolerance
                if(fabs(fi) < tolerance){
                    data_out<<"root: x = "<<xi<<endl;
                    break;
                }
            }
        }
    }
	// Part (ii)
    data_out<<endl<<"Roots found for part (ii) "<<endl;
    for(step = -5; step < 5; step+=0.5){
        x_left = step;
        x_right = step+0.5;
        f_left = f2(x_left);
        f_right = f2(x_right);
        if(f_left * f_right < 0) {          // opposite signs
            xi = (x_left + x_right)/2;		// guess for Newton root finding: average between left and right brackets
            // Loop through Newton-Raphson
            for(i = 0; i < N; i++){
                x = xi;
                f = f2(x);
                df = df2(x);
                // Newton step
                xi = x - f/df;
                fi = f2(xi);
                // Check if convergence is within tolerance
                if(fabs(fi) < tolerance){
                    data_out<<"root: x = "<<xi<<endl;
                    break;
                }
            }
        }
    }
    return 0;
}