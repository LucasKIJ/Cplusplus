#ifndef MATTESTSDEF
#define MATTESTSDEF
#include "Matrix.hpp"
namespace MatrixTests
{
    void RunTests();
    void Equality();
    void Add();
    void Minus();
    void MultiplySquare();
    void MultiplyRect();
    void MultiplyConst();
    void UnaryMinus();
    void Transpose();
    void RowMax();
    void ColMax();
    void SwapRow();
    void SwapCol();
    void Eye();
    void Swap();
    void LU();
    void Det();
    void QR();
    void HessReduction();
    void EigenVal();
    
    void NormAll();
    void Size();
    void GaussianElimination();
}
#endif