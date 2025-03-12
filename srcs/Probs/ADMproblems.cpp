#include <Geodesics.h>
extern double (*geodesic_points)[5];
extern int num_points;
extern double a;

int grid_setup() {
    Grid grid_obj;
	
	grid_obj.allocateGlobalGrid();
	grid_obj.initializeKerrData(grid_obj);
	grid_obj.evolve(grid_obj, 0.0000001, 4);
	printf("end of compute\n");
    return 0;
}
