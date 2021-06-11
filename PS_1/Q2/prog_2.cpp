#include <iostream>
#include <string>
int main()
{
    int n, x, prod;
    prod = 1;
    std::cout<<"How many integers do you want to multiply \n";
    std::cin>>n;
    for (int i =0; i < n; i++) {
        std::cout<<"Enter integer \n";
        std::cin>>x;
        prod =  prod*x;
    }
    std::cout<<"The product is: \n";
    std::cout<<prod;
}