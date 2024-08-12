#include <stdio.h>
#include <stdlib.h>

void rk(unsigned n, void fun(unsigned, double, double*, double*), 
        double h, double x, double* y, double* f, double* s, double* yk)
{
int i;
double xk;
//double f[N_EQ], s[N_EQ], yk[N_EQ];

fun(n, x, y, f);

for(i = 0; i < n; ++i)
   { s[i] = f[i]; yk[i] = y[i] + (h/2.)*f[i]; }
xk = x + h/2.;

fun(n, xk, yk, f);

for(i = 0; i < n; ++i) 
   { s[i] = s[i] + 2.*f[i]; yk[i] = y[i] + (h/2.)*f[i]; }
fun(n, xk, yk, f);

for(i = 0; i < n; ++i) 

   { s[i] = s[i] + 2.*f[i]; yk[i] = y[i] + h*f[i]; }
xk = x + h;
fun(n, xk, yk, f);

for(i = 0; i < n; ++i) 
   { y[i] = y[i] + (h/6.)*(s[i] + f[i]); }

}
