#include <iostream>
#include <cmath>
using namespace std;

// Define the function to solve
double func(double x, double y, double z) {
    return -((x * x - 2) * y) / (x * x);
}

// Define the derivative of the function
double deriv(double x, double y, double z) {
    return z;
}

// Define the 4th order Runge-Kutta method
void rk4(double x0, double y0, double z0, double xf, double h) {
    double x = x0, y = y0, z = z0;
    cout << "RungeKut method" << "\n";
    while (x < xf) {
        // Calculate the k values
        double k1y = h * deriv(x, y, z);
        double k1z = h * func(x, y, z);
        double k2y = h * deriv(x + h / 2, y + k1y / 2, z + k1z / 2);
        double k2z = h * func(x + h / 2, y + k1y / 2, z + k1z / 2);
        double k3y = h * deriv(x + h / 2, y + k2y / 2, z + k2z / 2);
        double k3z = h * func(x + h / 2, y + k2y / 2, z + k2z / 2);
        double k4y = h * deriv(x + h, y + k3y, z + k3z);
        double k4z = h * func(x + h, y + k3y, z + k3z);

        // Calculate the new values of y and z
        y = y + (k1y + 2 * k2y + 2 * k3y + k4y) / 6;
        z = z + (k1z + 2 * k2z + 2 * k3z + k4z) / 6;

        // Update the value of x
        x = x + h;

        // Output the current values of x and y
        cout << "x: " << x << " y: " << y << endl;
    }
}

void euler(double x0, double y0, double z0, double xf, double h) {
    double x = x0, y = y0, z = z0;
    cout << "Euler's method" << "\n";
    while (x < xf) {
        // Calculate the intermediate values
        double y1 = y + h * deriv(x, y, z);
        double z1 = z + h * func(x, y, z);
        double y2 = y + h * deriv(x + h / 2, y1, z1);
        double z2 = z + h * func(x + h / 2, y1, z1);
        double y3 = y + h * deriv(x + h / 2, y2, z2);
        double z3 = z + h * func(x + h / 2, y2, z2);
        double y4 = y + h * deriv(x + h, y3, z3);
        double z4 = z + h * func(x + h, y3, z3);

        // Calculate the new values of y and z
        y = y + (h / 6) * (deriv(x, y, z) + 2 * deriv(x + h / 2, y1, z1) + 2 * deriv(x + h / 2, y2, z2) + deriv(x + h, y3, z3));
        z = z + (h / 6) * (func(x, y, z) + 2 * func(x + h / 2, y1, z1) + 2 * func(x + h / 2, y2, z2) + func(x + h, y3, z3));

        // Update the value of x
        x = x + h;

        // Output the current values of x and y
        cout << "x: " << x << " y: " << y << endl;
    }
}



// Main function to run the Runge-Kutta method
int main() {
    double x0 = 1.0, y0 = 1.0, z0 = 0.0, xf = 2.0, h = 0.1;

    euler(x0, y0, z0, xf, h);
    rk4(x0, y0, z0, xf, h);

    return 0;
}
