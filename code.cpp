#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;

// ---------- Function for Lagrange Interpolation ----------
double lagrangeInterpolation(vector<double> x, vector<double> y, double xp) {
    double yp = 0;
    int n = x.size();
    for (int i = 0; i < n; i++) {
        double p = 1;
        for (int j = 0; j < n; j++) {
            if (i != j)
                p *= (xp - x[j]) / (x[i] - x[j]);
        }
        yp += p * y[i];
    }
    return yp;
}

// ---------- Function for Newton’s Forward Interpolation ----------
double newtonForward(vector<double> x, vector<double> y, double xp) {
    int n = x.size();
    double h = x[1] - x[0];

    // Create forward difference table
    double diff[10][10];
    for (int i = 0; i < n; i++)
        diff[i][0] = y[i];

    for (int j = 1; j < n; j++) {
        for (int i = 0; i < n - j; i++) {
            diff[i][j] = diff[i + 1][j - 1] - diff[i][j - 1];
        }
    }

    // Apply Newton forward formula
    double s = (xp - x[0]) / h;
    double yp = y[0];
    double term = 1;

    for (int i = 1; i < n; i++) {
        term *= (s - (i - 1)) / i;
        yp += term * diff[0][i];
    }

    return yp;
}

// ---------- Function for Newton’s Backward Interpolation ----------
double newtonBackward(vector<double> x, vector<double> y, double xp) {
    int n = x.size();
    double h = x[1] - x[0];

    // Create backward difference table
    double diff[10][10];
    for (int i = 0; i < n; i++)
        diff[i][0] = y[i];

    for (int j = 1; j < n; j++) {
        for (int i = n - 1; i >= j; i--) {
            diff[i][j] = diff[i][j - 1] - diff[i - 1][j - 1];
        }
    }

    // Apply Newton backward formula
    double s = (xp - x[n - 1]) / h;
    double yp = y[n - 1];
    double term = 1;

    for (int i = 1; i < n; i++) {
        term *= (s + (i - 1)) / i;
        yp += term * diff[n - 1][i];
    }

    return yp;
}

int main() {
    cout << fixed << setprecision(4);

    // ---------- Q1: Lagrange’s Interpolation ----------
    vector<double> x1 = {5, 7, 11, 13, 17};
    vector<double> y1 = {150, 392, 1452, 2366, 5202};
    double x_val = 10;

    double y_val = lagrangeInterpolation(x1, y1, x_val);
    cout << "Q1: Using Lagrange’s Interpolation\n";
    cout << "Value of y when x = " << x_val << " is " << y_val << "\n\n";

    // ---------- Q2: Newton’s Interpolation ----------
    vector<double> x2 = {1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3};
    vector<double> y2 = {5.474, 6.050, 6.686, 7.389, 8.166, 9.025, 9.974};

    double x_forward = 1.85;
    double x_backward = 2.25;

    double y_forward = newtonForward(x2, y2, x_forward);
    double y_backward = newtonBackward(x2, y2, x_backward);

    cout << "Q2: Using Newton’s Interpolation\n";
    cout << "Value of y when x = " << x_forward << " (Forward) is " << y_forward << "\n";
    cout << "Value of y when x = " << x_backward << " (Backward) is " << y_backward << "\n";

    return 0;
}
