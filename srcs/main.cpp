#include <Geodesics.h>

double (*geodesic_points)[5] = NULL;
int num_points = 0;
double a = 0.935;

int main(int argc, char **argv)
{
	if (argc < 4) {
		printf("Usage: <options>\n");
		printf("Options:\n");
		printf("       -G <Spin value a> - Geodesic calculation\n");
		printf("       -L <Spin value a> - Light geodesics calculation\n");
		printf("       -R <Spin value a> - Riemann tensor calculation\n");
		printf("       -M <Spin value a> - Metric tensor calculation (a = 0 -> Schwarzschild or a > 0 -> Kerr)\n");
		printf("       -S <Spin value a> - Black hole shadow generation\n");	
		printf("       -C <Spin value a> - ADM solver Kerr-Schild coordinates\n");
		return 0;
	}

	a = atof(argv[2]);
	if (strcmp(argv[1], "-R") == 0) {
		if (argc < 4) {
			printf("Usage: -R <Spin value a> <Metric>\n");
			printf("Metric: schwarzschild, kerr, kerr-newman, kds\n");
			return 0;
		} 
		Riemann_tensor(argv[3]);
	} else if (strcmp(argv[1], "-G") == 0) {
		if (argc < 3) {
			printf("Usage: -G <Spin value a>\n");
			return 0;
		} 
		Geodesics_prob();	
	}else if (strncmp(argv[1], "-L", 2) == 0) {
		if (argc < 3) {
			printf("Usage: -L <Spin value a>\n");
			return 0;
		}
		light_geodesics_prob();
	} else if (strncmp(argv[1], "-M", 2) == 0) {
		if (argc < 3) {
			printf("Usage: -M <Spin value a>\n");
			return 0;
		}
		Metric_prob();
	} else if (strncmp(argv[1], "-C", 2) == 0) {
		if (argc < 3) {
			printf("Usage: -C <Spin value a>\n");
			printf("ADM solver Kerr-Schild coordinates\n");
			printf("Metric: schwarzschild, kerr or Minkowski\n");
			return 0;
		}
		grid_setup();
	}else {
		printf("Invalid option\n");
		return 0;
	}
}
