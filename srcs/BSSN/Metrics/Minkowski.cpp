#include <Geodesics.h>

void Grid::initializeData_Minkowski()
{
    printf("\n=== Initialisation Minkowski ===\n");

    double x_min = -128.0, x_max = 128.0;
    double y_min = -128.0, y_max = 128.0;
    double z_min = -128.0, z_max = 128.0;
    double dx = (x_max - x_min)/(NX-1);
    double dy = (y_max - y_min)/(NY-1);
    double dz = (z_max - z_min)/(NZ-1);

    for(int i=0; i<NX; i++)
    {
        double x = x_min + i*dx;
        for(int j=0; j<NY; j++)
        {
            double y = y_min + j*dy;
            for(int k=0; k<NZ; k++)
            {
                Cell2D &cell = globalGrid[i][j][k];

                cell.gauge.alpha = 1.0;
                cell.gauge.beta[0] = 0.0;
                cell.gauge.beta[1] = 0.0;
                cell.gauge.beta[2] = 0.0;

                for(int a=0; a<3; a++)
                {
                    for(int b=0; b<3; b++)
                    {
                        cell.geom.gamma[a][b] = (a == b) ? 1.0 : 0.0;
                    }
                }

                for(int a=0; a<3; a++){
                    for(int b=0; b<3; b++){
                        cell.curv.K[a][b] = 0.0;
                    }
                }
            }
        }
    }

    printf("Vérification Minkowski sur quelques points:\n");
    for(int test_i=0; test_i<3; test_i++)
    {
        for(int test_j=0; test_j<3; test_j++)
        {
            for(int test_k=0; test_k<3; test_k++)
            {
                Cell2D &cell = globalGrid[test_i][test_j][test_k];
                printf("Point (%d,%d,%d): alpha=%f, geom.gamma=\n", test_i, test_j, test_k, cell.gauge.alpha);
                for(int a=0; a<3; a++)
                {
                    printf("  ");
                    for(int b=0; b<3; b++)
                    {
                        printf("%f ", cell.geom.gamma[a][b]);
                    }
                    printf("\n");
                }
            }
        }
    }
    printf("=== Fin initialisation Minkowski ===\n");
}

