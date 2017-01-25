// LU.cpp
//
// This program decomposes a matrix into LU format and solves a system of equations given a b vector
// 
//*****************
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <stdio.h>
using namespace std;

// the main part of the program - should always be type "int"

int main() {

	// Variable declarations
    int i,j,k; // loop counters
	const int degree = 4; // degree of polynomial fit
	const int n = degree+1; // size of matrices
	const int m = 22; // number of data points
	double sum; // dummy variable to keep track of sums in decomposition

	// Our x-axis values
	double x_T[] = {25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 150, 200, 250, 273, 300, 350, 400, 450, 500, 600, 700};

	// Our y-axis values
	double y_k[] = {459.5, 317.7, 233.3, 180.3, 143.2, 118.4, 85.2, 66.6, 54.4, 46.1, 40.2, 23.3, 16.5, 12.6, 11.5, 10.5, 8.8, 7.6, 6.9, 6.1, 5.1, 4.5};

	// Solve Ax=b using LU decomposition
	double A[n][n]; // coefficient matrix
    double L[n][n]; // lower triangular matrix
    double U[n][n]; // upper triangular matrix
	double x[n]; // solution vector
    double y[n]; // intermediate solution vector
	double b[n];// right side vector

	//		create array of x summations
	double sum_x[2*degree]; // array of x summations
	
	for (k=0; k<2*degree; k++) {
		sum_x[k] = 0; // intial value of xi sum
		for (i=0; i<m; i++) {
			sum_x[k] += pow(x_T[i], k+1);
		}
	}

	//		build the b vector
	for (k=0; k<n; k++) {
		b[k] = 0; // intial value of (xi^n)*yi sum
		for (i=0; i<m; i++) {
			b[k] += pow(x_T[i], k)*y_k[i];
		}
		cout<<b[k]<<'\n';  // print b vector values
	}

	//		build A matrix
	A[0][0] = m;
	for (k=0; k<2*degree; k++) {
		for (i=0; i<n; i++) {
			for (j=0; j<n; j++) {
				if (i+j==k+1) {
					A[i][j]=sum_x[k];
				}
				cout<<A[i][j]<<"\t";
			}
		}
	}

	//		print out the A Matrix
	for(i=0;i<n;i++) {  // loop n times for the n rows
        for(j=0;j<n;j++)  // loop n times for the n cols
        {
            cout<<A[i][j]<<"\t";  // display the current element out of the array
        }
    cout<<endl;  // when the inner loop is done, go to a new line
	}

	// Copy Pasta LU Code
	// initialize L and U
   
   for(i=0;i<n;i++) {
	   for(j=0;j<n;j++){
		   L[i][j] = 0.;
		   U[i][j] = 0.;
	   }
   }
    
   for(i=0;i<n;i++) U[i][i] = 1.; // for Crout reduction
  
   // Do the LU decomposition using the Crout reduction
   
   for(i=0;i<n;i++) { // loop over pairs of L columns and U rows. The three levels of loop indicate n^3 behavior
	  
     // first, column i of L
     
   		for(j=i;j<n;j++) { // i is the row
			sum=0.;
			for(k=0;k<=i-1;k++) sum = sum + L[j][k]*U[k][i];
			L[j][i] = A[j][i] - sum;
		}
		
	// second, row i of U
		   
		  for(j=i+1;j<n;j++){ // j is the column
			   sum = 0.;
			   for(k=0;k<=i-1;k++) sum = sum + L[i][k]*U[k][j];
			   U[i][j] = (A[i][j] - sum)/L[i][i];
		   }
   }
   
   // output intermediate data to screen
   
   for(i=0;i<n;i++) {
	   for(j=0;j<n;j++) cout<<A[i][j]<<'\t';
	   cout<<endl;
   }
   cout<<endl;		   
    
   for(i=0;i<n;i++) {
	   for(j=0;j<n;j++) cout<<L[i][j]<<'\t';
	   cout<<endl;
   }
   cout<<endl;	
   
   for(i=0;i<n;i++) {
	   for(j=0;j<n;j++) cout<<U[i][j]<<'\t';
	   cout<<endl;
   }
   cout<<endl;	
   
   // solve the system of equations
   
   // could loop over a series of b vectors here once the decomposition has been done
   
   // first, find the y vector
   
   y[0] = b[0]/L[0][0];
   for(i=1;i<n;i++) {
	   sum = 0.;
	   for(j=0;j<i;j++) sum = sum + L[i][j]*y[j];
	   y[i] = (b[i]-sum)/L[i][i];
   }
   
   // second, find the x vector
   
   x[n-1] = y[n-1];
   for(i=1;i<n;i++) {
	   j = n-i-1;
	   sum = 0;
	   for(k=j+1;k<n;k++) sum = sum + U[j][k]*x[k];
	   x[j] = y[j] - sum;
   }
   
   // output data
   for(i=0;i<n;i++) cout<<x[i]<<'\t';
   cout<<endl;
}