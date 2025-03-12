#pragma once

#include <Grid.h>

class Grid;

class Plot{
	public:
		Grid* grid;
		friend class Grid;
		Plot() = default;
		~Plot() = default;

		void export_gamma_slice(int j);
		void export_K_slice(int j);
		void export_alpha_slice(int j);
		void export_gauge_slice(int j);
};
