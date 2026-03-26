#pragma once

#include <string>

class Grid;

void logger_evolve(const Grid &grid_obj, float dt, int nstep);
void appendConstraintL2ToCSV(const Grid &grid_obj, const std::string &filename,
                             float time);
void export_Atildedt_slide(Grid &grid_obj, float time);
void export_K_3D(Grid &grid_obj);
void export_alpha_slice(Grid &grid_obj, int j);
void export_gauge_slice(Grid &grid_obj, int j);
void export_K_slice(Grid &grid_obj, int j);
void export_gamma_slice(Grid &grid_obj, int j, float time);
void export_tilde_gamma_3D(Grid &grid_obj);
void export_gauge_slice_anim(Grid &grid_obj, int j, float time);
void export_gauge_slice_xy(Grid &grid_obj, int k, float time);
