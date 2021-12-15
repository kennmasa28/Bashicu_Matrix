#include "athena_arrays.hpp"

class Matrix{
    public:
    int WhetherZero(AthenaArray<int>&, int, int);
    int ActivateFunction(int,int);
    int NonZeroLowLine(AthenaArray<int>&, int, int);
    int DecideBadRoot(AthenaArray<int>&, int, int);
    int IsColumnOfDirectAncestor(AthenaArray<int>&, int, int, int);
    int ColumnOfParent(AthenaArray<int>&, int, int);
    void DecideGoodPart(AthenaArray<int>&, AthenaArray<int>&, int, int);
    void DecideDelta(AthenaArray<int>&, AthenaArray<int>&, int, int, int, int);
    void DecideA(AthenaArray<int>&, AthenaArray<int>&, int, int, int, int);
    void DecideBadPart(AthenaArray<int>&, AthenaArray<int>&,  AthenaArray<int>&, AthenaArray<int>&, int, int, int, int);
    void NewMatrix(AthenaArray<int>&,  AthenaArray<int>&, AthenaArray<int>&, int, int, int, int);
};

class InOut{
    public:
    void GetLineAndColumn(AthenaArray<int>&);
    void InputMatrix(AthenaArray<int>&);
    std::vector<std::string> split(std::string, std::string);
    void OutputMatrix(AthenaArray<int>&, int, int,int);
    void Outputb(AthenaArray<int>&, int, int, int);
};

