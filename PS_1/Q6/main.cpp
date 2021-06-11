#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>

int main()
{
    double x, y, max_err, err;
    std::ifstream input("../Q1_5/xy.dat");
    assert(input.is_open());
    max_err = 0;
    while (true){
        input >> x >> y;
        if (input.eof()){
            break;
        }
        err = abs(exp(-x) - y);
        if (max_err < err){
            max_err = err;
        }
    }
    std::cout << "The maximum absolute error is: " << max_err;
}