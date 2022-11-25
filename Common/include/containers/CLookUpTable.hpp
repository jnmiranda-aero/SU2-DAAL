/*!
 * \file CLookupTable.hpp
 * \brief tabulation of fluid properties
 * \author D. Mayer, T. Economon, E.C.Bunschoten
 * \version 7.4.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <iomanip>
#include <string>
#include <vector>
#include <array>

#include "../../Common/include/option_structure.hpp"
#include "CFileReaderLUT.hpp"
#include "CTrapezoidalMap.hpp"

/*!
 * \brief Look up table.
 * \ingroup LookUpInterp
 */
class CLookUpTable {
 protected:
  int rank; /*!< \brief MPI Rank. */

  std::string file_name_lut;
  std::string version_lut;
  std::string version_reader;
  std::string name_CV1,   // Name of controlling variable 1 
              name_CV2;   // Name of xontrolling variable 2

  unsigned long *n_points;
  unsigned long *n_triangles;
  unsigned long n_variables;
  unsigned long *n_hull_points;
  unsigned long n_table_levels = 1;
  su2double *z_values_levels = nullptr;

  /*! 
   * \brief the lower and upper limits of the z, y and x variable.
   */
  std::pair<su2double, su2double> limits_table_z;
  std::pair<su2double, su2double> *limits_table_y;
  std::pair<su2double, su2double> *limits_table_x;

  /*! \brief Holds the variable names stored in the table file.
   * Order is in sync with data
   */
  su2vector<std::string> names_var;

  /*! \brief
   * Holds all data stored in the table. First index addresses the variable
   * while second index addresses the point.
   */
  su2activematrix *table_data;

  su2matrix<unsigned long> *triangles;

  std::vector<su2vector<unsigned long> > *edges;
  su2vector<std::vector<unsigned long> > *edge_to_triangle;

  /*! \brief 
   * The hull contains the boundary of the lookup table.
   */
  su2vector<unsigned long> *hull;

  CTrapezoidalMap *trap_map_x_y;

  /*! \brief 
   * vector of all the weight factors for the interpolation.
   */
  su2vector<su2activematrix> *interp_mat_inv_x_y;

  /*! \brief 
   * returns the index to the variable in the lookup table.
   */
  inline unsigned int GetIndexOfVar(const std::string& nameVar) const {
    int index = find(names_var.begin(), names_var.end(), nameVar) - names_var.begin();
    if (index == int(names_var.size())) {
      index = -1;
      std::string error_msg = "Variable '";
      error_msg.append(nameVar);
      error_msg.append("' is not in the lookup table.");
      SU2_MPI::Error(error_msg, CURRENT_FUNCTION);
    }
    return index;
  }

  /*!
   * \brief Get the pointer to the column data of the table (density, temperature, source terms, ...).
   * \returns pointer to the column data.
   */
  inline const su2double* GetDataP(const std::string& name_var, unsigned long i_level=0) const {
    return table_data[i_level][GetIndexOfVar(name_var)];
  }

  /*!
   * \brief find the table limits, i.e. the minimum and maximum values of the 2 independent.
   * controlling variables. We put the values in the variables.
   * limits_table_x[2] and limit_table_y[2].
   * \param[in] name_x - the string name for the first controlling variable.
   * \param[in] name_y - the string name of the second controlling variable.
   */
  void FindTableLimits(const std::string& name_x, const std::string& name_y);

  /*!
   * \brief construct a list of all the edges and a list of the pair of elements left and right of the edge.
   */
  void IdentifyUniqueEdges();

  /*!
   * \brief read the lookup table from file and store the data.
   * \param[in] file_name_lut - the filename of the lookup table.
  */
  void LoadTableRaw(const std::string& file_name_lut);

  /*!
   * \brief compute vector of all (inverse) interpolation coefficients "interp_mat_inv_x_y" of all triangles.
  */
  void ComputeInterpCoeffs();

  /*!
   * \brief compute the inverse matrix for interpolation.
   * \param[in] vec_x - pointer to first coordinate (progress variable).
   * \param[in] vec_y - pointer to second coordinate (enthalpy).
   * \param[in] point_ids - single triangle data.
   * \param[out] interp_mat_inv - inverse matrix for interpolation.
  */
  void GetInterpMatInv(const su2double* vec_x, const su2double* vec_y, std::array<unsigned long,3>& point_ids,
                       su2activematrix& interp_mat_inv);

  /*!
   * \brief compute the interpolation coefficients for the triangular interpolation.
   * \param[in] val_x - value of first coordinate (progress variable).
   * \param[in] val_y - value of second coordinate (enthalpy).
   * \param[in] interp_mat_inv - inverse matrix for interpolation.
   * \param[out] interp_coeffs - interpolation coefficients.
  */
  void GetInterpCoeffs(su2double val_x, su2double val_y, su2activematrix& interp_mat_inv,
                       std::array<su2double,3>& interp_coeffs);

  /*!
   * \brief compute interpolated value of a point P in the triangle.
   * \param[in] val_samples - pointer to the variable data.
   * \param[in] val_triangle - ID to the triangle.
   * \param[in] val_interp_coeffs - interpolation coefficients using the point data in P.
   * \returns resulting value of the interpolation.
  */
  su2double Interpolate(const su2double* val_samples, std::array<unsigned long,3>& val_triangle,
                        std::array<su2double,3>& val_interp_coeffs);

  void Linear_Interpolation(su2double val_z, std::pair<unsigned long, unsigned long> &inclusion_levels, su2double&lower_value,su2double&upper_value, su2double*&var_vals);


  void Linear_Interpolation(su2double val_z, std::pair<unsigned long, unsigned long> &inclusion_levels, std::vector<su2double>&lower_values,std::vector<su2double>&upper_values, std::vector<su2double*>&var_vals);

  /*!
   * \brief find the point on the hull (boundary of the table) that is closest to the point P(val_x,val_y).
   * \param[in] val_x - first coordinate of point P(val_x,val_y) to check.
   * \param[in] val_y - second coordinate of point P(val_x,val_y) to check.
   * \returns point id of the nearest neighbor on the hull.
   */
  unsigned long FindNearestNeighborOnHull(su2double val_x, su2double val_y, unsigned long i_level=0);

  /*!
   * \brief determine if a point P(val_x,val_y) is inside the triangle val_id_triangle.
   * \param[in] val_x - first coordinate of point P(val_x,val_y) to check.
   * \param[in] val_y - second coordinate of point P(val_x,val_y) to check.
   * \param[in] val_id_triangle - ID of the triangle to check.
   * \returns true if the point is in the triangle, false if it is outside.
   */
  bool IsInTriangle(su2double val_x, su2double val_y, unsigned long val_id_triangle, unsigned long i_level=0);

  /*!
   * \brief compute the area of a triangle given the 3 points of the triangle.
   * \param[in] x1 - the coordinates of the points P1(x1,y1), P2(x2,y2) and P3(x3,y3).
   * \param[in] y1 - the coordinates of the points P1(x1,y1), P2(x2,y2) and P3(x3,y3).
   * \param[in] x2 - the coordinates of the points P1(x1,y1), P2(x2,y2) and P3(x3,y3).
   * \param[in] y2 - the coordinates of the points P1(x1,y1), P2(x2,y2) and P3(x3,y3).
   * \param[in] x3 - the coordinates of the points P1(x1,y1), P2(x2,y2) and P3(x3,y3).
   * \param[in] y3 - the coordinates of the points P1(x1,y1), P2(x2,y2) and P3(x3,y3).
   * \returns the absolute value of the area of the triangle.
   */
  inline su2double TriArea(su2double x1, su2double y1, su2double x2, su2double y2, su2double x3, su2double y3) {
    return abs((x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) * 0.5);
  }

 public:
  CLookUpTable(const std::string& file_name_lut, const std::string& name_x, const std::string& name_y);
  ~CLookUpTable();

  /*!
   * \brief print information to screen.
   */
  void PrintTableInfo();

  /*!
   * \brief lookup 1 value of the single variable "val_name_var" using controlling variable values(val_x,val_y).
   * \param[in] val_name_var - string name of the variable to look up.
   * \param[out] val_var - the stored value of the variable to look up.
   * \param[in] val_x - value of controlling variable 1.
   * \param[in] val_y - value of controlling variable 2.
   * \returns 1 if the lookup and subsequent interpolation was a success, 0 if not.
   */
  unsigned long LookUp_XY(const std::string& val_name_var, su2double* val_var, su2double val_x, su2double val_y, unsigned long i_level=0);

  /*!
   * \brief lookup 1 value for each of the variables in "val_name_var" using controlling variable values(val_x,val_y).
   * \param[in] val_names_var - vector of string names of the variables to look up.
   * \param[out] val_vars - pointer to the vector of stored values of the variables to look up.
   * \param[in] val_x - value of controlling variable 1.
   * \param[in] val_y - value of controlling variable 2.
   * \returns 1 if the lookup and subsequent interpolation was a success, 0 if not.
   */
  unsigned long LookUp_XY(const std::vector<std::string>& val_names_var, std::vector<su2double*>& val_vars, su2double val_x,
                                su2double val_y, unsigned long i_level=0);

  /*!
   * \brief lookup the value of the variable "val_name_var" using controlling variable values(val_x,val_y).
   * \param[in] val_name_var - string name of the variable to look up.
   * \param[out] val_var - the stored value of the variable to look up.
   * \param[in] val_x - value of controlling variable 1.
   * \param[in] val_y - value of controlling variable 2.
   * \returns 1 if the lookup and subsequent interpolation was a success, 0 if not.
   */
  unsigned long LookUp_XY(const std::vector<std::string>& val_names_var, std::vector<su2double>& val_vars, su2double val_x,
                                su2double val_y, unsigned long i_level=0);

  /*!
   * \brief lookup the value of the variable "val_name_var" using controlling variable values(val_x,val_y,val_z).
   * \param[in] val_name_var - string name of the variable to look up.
   * \param[out] val_var - the stored value of the variable to look up.
   * \param[in] val_x - value of controlling variable 1.
   * \param[in] val_y - value of controlling variable 2.
   * \param[in] val_z - value of controlling variable 3.
   * \returns 1 if the lookup and subsequent interpolation was a success, 0 if not.
   */
  unsigned long LookUp_XYZ(const std::string& val_name_var, su2double* val_var, su2double val_x, su2double val_y, su2double val_z);

  /*!
   * \brief lookup the value of the variable "val_name_var" using controlling variable values(val_x,val_y,val_z).
   * \param[in] val_name_var - string name of the variable to look up.
   * \param[out] val_var - the stored value of the variable to look up.
   * \param[in] val_x - value of controlling variable 1.
   * \param[in] val_y - value of controlling variable 2.
   * \param[in] val_z - value of controlling variable 3.
   * \returns 1 if the lookup and subsequent interpolation was a success, 0 if not.
   */
  unsigned long LookUp_XYZ(const std::vector<std::string>& val_names_var, std::vector<su2double*>& val_vars, su2double val_x,
                                su2double val_y, su2double val_z=0);
  
  std::pair<unsigned long, unsigned long> FindInclusionLevels(su2double val_z, bool *within_limits);

  /*!
   * \brief determine the minimum and maximum value of the second controlling variable.
   * \returns pair of minimum and maximum value of controlling variable 2.
   */
  inline std::pair<su2double, su2double> GetTableLimitsY(unsigned long i_level = 0) const {
    return limits_table_y[i_level];
  }

  /*!
   * \brief determine the minimum and maximum value of the first controlling variable.
   * \returns pair of minimum and maximum value of controlling variable 1.
   */
  inline std::pair<su2double, su2double> GetTableLimitsX(unsigned long i_level = 0) const {
    return limits_table_x[i_level];
  }

  
};