#include <iostream>
#include <fstream>
#include <cassert>
int main()
{
    std::ifstream input("numbers.dat");
    assert(input.is_open()); //Edit 23/04/2018
    double x, y;
    while (true)
    {
        input >> x >> y;
        if (input.eof())
        {
            break;
        }
        std::cout << x << " " << y << "\n";
    }
    return 0;
}
