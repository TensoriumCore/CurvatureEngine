#pragma once
#include <Grid.h> 
class Grid;
class Derivatives
{
    public:
        Derivatives() = default;
        ~Derivatives() = default;
        friend class Grid;
        friend class GridTensor;
        Grid* grid;
        template <class T>
        inline double fourth_order_diff(const T &plus2, const T &plus1, 
                const T &minus1, const T &minus2, double dx)
        {
            return (-plus2 + 8.0*plus1 - 8.0*minus1 + minus2) / (12.0*dx);
        }
        template <class T>
        inline double second_order_diff(const T &plus1, const T &minus1, double dx)
        {
            return (plus1 - minus1)/(2.0*dx);
        }
        double partialX_alpha(int i, int j, int k);
        double partialY_alpha(int i, int j, int k);
        double partialZ_alpha(int i, int j, int k);

        double partialXX_alpha(int i, int j, int k);
        double partialYY_alpha(int i, int j, int k);
        double partialZZ_alpha(int i, int j, int k);

        double partialX_betacomp(int i, int j, int k, int comp);
        double partialY_betacomp(int i, int j, int k, int comp);
        double partialZ_betacomp(int i, int j, int k, int comp);

        double partialXY_alpha(int i, int j, int k);
        double partialXZ_alpha(int i, int j, int k);
        double partialYZ_alpha(int i, int j, int k);

        double second_partial_alpha(int i, int j, int k, int a, int b);
};
