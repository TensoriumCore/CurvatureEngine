#include <Geodesics.h>

int grid_setup() {
    Grid grid_obj;
	
	grid_obj.allocateGlobalGrid();
	grid_obj.initializeKerrData(grid_obj);
	grid_obj.evolve(grid_obj, 0.0001, 3);
	printf("end of compute\n");
    return 0;
}
