#include <iostream>
#include <string>
int main()
{
    int x, prod;
    prod = 1;
    while (true) {
        std::cout<<"Enter integer \n";
        std::cin>>x;
        if (x == -1){
            break;
        }
        prod =  prod*x;
    }
    std::cout<<"The product is: \n";
    std::cout<<prod;
}