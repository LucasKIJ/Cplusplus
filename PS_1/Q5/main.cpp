#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>

int main()
{
    int N;
    std::cout << "How many grid points\? \n";
    std::cin >> N;
    assert(N > 1);
    double y[N];
    double x[N];
    const double h = 1.0 / N;
    x[0] = 0;
    y[0] = 1;

    std::ofstream out("xy.dat");
    assert(out.is_open());
    out.setf(std::ios::scientific|std::ios::showpos);
    out << x[0] << " " << y[0] << "\n";

    for (int i = 0; i < N; i++){
        x[i+1] = x[i] + h;
        y[i+1] = y[i] / (1+h);
        out << x[i+1] << " " << y[i+1] << "\n";
    }
    out.close();
    return 0;
}