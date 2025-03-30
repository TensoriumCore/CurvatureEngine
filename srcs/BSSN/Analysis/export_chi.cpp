#include <Geodesics.h>

void Grid::export_BH_positions(double time, const std::string& filename) {
    struct BH {
        double chi;
        int i, j, k;
    };

    BH bh1 = {1e10, -1, -1, -1};
    BH bh2 = {1e10, -1, -1, -1};

    for (int i = 2; i < NX - 2; ++i) {
        for (int j = 2; j < NY - 2; ++j) {
            for (int k = 2; k < NZ - 2; ++k) {
                double chi_val = globalGrid[i][j][k].dt_chi;

                if (chi_val < bh1.chi) {
                    double dist = std::sqrt(
                        std::pow(i - bh1.i, 2) + std::pow(j - bh1.j, 2) + std::pow(k - bh1.k, 2));
                    if (dist > 5 || bh1.i == -1) {
                        bh2 = bh1;
                        bh1 = {chi_val, i, j, k};
                    }
                } else if (chi_val < bh2.chi) {
                    double dist = std::sqrt(
                        std::pow(i - bh1.i, 2) + std::pow(j - bh1.j, 2) + std::pow(k - bh1.k, 2));
                    if (dist > 5) {
                        bh2 = {chi_val, i, j, k};
                    }
                }
            }
        }
    }

    std::ofstream out(filename, std::ios::app); 
    if (!out.is_open()) {
        std::cerr << "Failed to open " << filename << " for writing.\n";
        return;
    }

    double x1 = -9 + bh1.i * DX;
    double y1 = -9 + bh1.j * DY;
    double z1 = -9 + bh1.k * DZ;

    double x2 = -9 + bh2.i * DX;
    double y2 = -9 + bh2.j * DY;
    double z2 = -9 + bh2.k * DZ;
	printf("BH1: %f %f %f\n", x1, y1, z1);
	printf("BH2: %f %f %f\n", x2, y2, z2);
    out << std::fixed << std::setprecision(8)
        << time << ","
        << x1 << "," << y1 << "," << z1 << ","
        << x2 << "," << y2 << "," << z2 << "\n";

    out.close();
}

