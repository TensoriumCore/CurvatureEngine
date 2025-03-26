#include <Geodesics.h>

int grid_setup() {
    Grid grid_obj;
	
	grid_obj.allocateGlobalGrid();
	grid_obj.initializeKerrData(grid_obj);
	grid_obj.evolve(grid_obj, 0.0000001, 10);
	printf("end of compute\n");
    return 0;
}
