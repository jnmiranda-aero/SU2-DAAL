
/*!
 * \file CfluidFlamelet.cpp
 * \brief Main subroutines of CFluidFlamelet class
 * \author D. Mayer, T. Economon, N. Beishuizen
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

#include "../include/fluid/CFluidFlamelet.hpp"
#include "../../../Common/include/containers/CLookUpTable.hpp"

CFluidFlamelet::CFluidFlamelet(CConfig* config, su2double value_pressure_operating) : CFluidModel() {
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /* -- number of auxiliary species transport equations: 1=CO, 2=NOx --- */
  n_reactants = config->GetNReactants();
  n_CV = 2;
  n_scalars = n_CV + n_reactants;

  if (rank == MASTER_NODE) {
    cout << "n_scalars = " << n_scalars << endl;
    cout << "n_CV = " << n_CV << endl;
    cout << "n_reactants = " << n_reactants << endl;
  }

  config->SetNControllingVars(n_CV);

  config->SetNScalars(n_scalars);

  if (rank == MASTER_NODE) {
    cout << "*****************************************" << endl;
    cout << "***   initializing the lookup table   ***" << endl;
    cout << "*****************************************" << endl;
  }

  table_scalar_names.resize(n_scalars);
  table_scalar_names[I_ENTH] = "EnthalpyTot";
  table_scalar_names[I_PROGVAR] = "ProgressVariable";
  /*--- auxiliary species transport equations---*/
  for(size_t i_aux=0; i_aux<n_reactants; i_aux++){
    table_scalar_names[n_CV + i_aux] = "Y-"+config->GetReactantName(i_aux);
  }

  config->SetLUTScalarNames(table_scalar_names);

  /*--- we currently only need one source term from the LUT for the progress variable
        and one additional source term for the each of the auxiliary species transport equations ---*/
  n_table_sources = 1 + n_reactants;
  config->SetNLUTSources(n_table_sources);

  table_source_names.resize(n_table_sources);
  table_sources.resize(n_table_sources);
  table_source_names[I_SRC_TOT_PROGVAR] = "ProdRateTot-PV";

  /*--- No source term for enthalpy ---*/

  /*--- source terms for auxiliary species transport equations ---*/
  for(size_t i_aux=0; i_aux<n_reactants; i_aux++){
    table_source_names[1 + i_aux] = "ProdRateTot-"+config->GetReactantName(i_aux);
  }

  config->SetLUTSourceNames(table_source_names);

  look_up_table =
      new CLookUpTable(config->GetFileNameLUT(), table_scalar_names[I_PROGVAR], table_scalar_names[I_ENTH]);

  n_lookups = config->GetNLookups();
  table_lookup_names.resize(n_lookups);
  for (int i_lookup = 0; i_lookup < n_lookups; ++i_lookup) {
    table_lookup_names.at(i_lookup) = config->GetLUTLookupName(i_lookup);
  }

  source_scalar.resize(n_scalars);
  lookup_scalar.resize(n_lookups);

  Pressure = value_pressure_operating;

  PreprocessLookUp();
}

CFluidFlamelet::~CFluidFlamelet() {
  if (look_up_table != NULL) delete look_up_table;
}

/* do a lookup for the list of variables in table_lookup_names, for visualization purposes */
unsigned long CFluidFlamelet::SetScalarLookups(su2double* val_scalars) {

  su2double enth = val_scalars[I_ENTH];
  su2double prog = val_scalars[I_PROGVAR];

  string name_enth = table_scalar_names[I_ENTH];
  string name_prog = table_scalar_names[I_PROGVAR];

  /* perform table look ups */
  unsigned long exit_code =
      look_up_table->LookUp_ProgEnth(table_lookup_names, lookup_scalar, prog, enth, name_prog, name_enth);

  return exit_code;
}

/* set the source terms for the transport equations*/
unsigned long CFluidFlamelet::SetScalarSources(su2double* val_scalars) {

  table_sources[0] = 0.0;

  string name_enth = table_scalar_names[I_ENTH];
  string name_prog = table_scalar_names[I_PROGVAR];

  /*--- value for the progress variable and enthalpy ---*/
  su2double enth = val_scalars[I_ENTH];
  su2double prog = val_scalars[I_PROGVAR];

  /* perform table look ups */
  unsigned long exit_code =
      look_up_table->LookUp_ProgEnth(varnames_Sources, val_vars_Sources, prog, enth, name_prog, name_enth);

  /*--- the source term for the progress variable is always positive, but we clip it just to be sure --- */
  source_scalar[I_PROGVAR] = max(EPS, table_sources[I_SRC_TOT_PROGVAR]);
  source_scalar[I_ENTH] = 0.0;
  /*--- source term for the auxiliary species transport equations---*/
  for(size_t i_aux = 0; i_aux < n_reactants; i_aux++) {
    source_scalar[n_CV + i_aux] = table_sources[1 + i_aux];
  }

  return exit_code;
}

void CFluidFlamelet::SetTDState_T(su2double val_temperature, const su2double* val_scalars) {
  su2double val_enth = val_scalars[I_ENTH];
  su2double val_prog = val_scalars[I_PROGVAR];

  string name_enth = table_scalar_names[I_ENTH];
  string name_prog = table_scalar_names[I_PROGVAR];

  /*--- add all quantities and their address to the look up vectors ---*/
  look_up_table->LookUp_ProgEnth(varnames_TD, val_vars_TD, val_prog, val_enth, name_prog, name_enth);

  /*--- compute Cv from Cp and molar weight of the mixture (ideal gas) ---*/
  Cv = Cp - UNIVERSAL_GAS_CONSTANT / molar_weight;
}

unsigned long CFluidFlamelet::GetEnthFromTemp(su2double* val_enth, su2double val_prog, su2double val_temp, su2double initial_value) {

  string name_prog = table_scalar_names[I_PROGVAR];
  string name_enth = table_scalar_names[I_ENTH];

  su2double   delta_temp_final = 0.01, /* convergence criterion for temperature in [K] */
              relaxation = 0.5,        /* Newton solver relaxation factor. */
              enth_iter = initial_value, /* Initial enthalpy value, default stetting is zero. */
              delta_enth,              /* Enthalpy residual. */
              delta_temp_iter;         /* Temperature residual. */

  unsigned long exit_code = 0,
                counter_limit = 50,
                counter = 0;

  bool converged = false;

  while (!converged && (counter++ < counter_limit)) {
    /* look up temperature and heat capacity */
    look_up_table->LookUp_ProgEnth(varnames_TD, val_vars_TD, val_prog, enth_iter, name_prog, name_enth);

    /* calculate delta_temperature */
    delta_temp_iter = val_temp - Temperature;
    if(abs(delta_temp_iter) < delta_temp_final){
      converged = true;
    }else{
      /* calculate delta_enthalpy following dh = cp * dT */
      delta_enth = Cp * delta_temp_iter;

      /* update enthalpy */
      enth_iter += relaxation * delta_enth;

      counter ++;
    }
  }

  /* set enthalpy value */
  *val_enth = enth_iter;

  if (!converged) {
    exit_code = 1;
  }

  return exit_code;
}

void CFluidFlamelet::PreprocessLookUp() {
  /*--- Set lookup names and variables for all relevant lookup processes in the fluid model ---*/

  /*--- Thermodynamic state variables ---*/
  varnames_TD.resize(7);
  val_vars_TD.resize(7);

  /* The string in varnames_TD is the actual string as it appears in the LUT file */
  varnames_TD[0] = "Temperature";
  val_vars_TD[0] = &Temperature;
  varnames_TD[1] = "Density";
  val_vars_TD[1] = &Density;
  varnames_TD[2] = "Cp";
  val_vars_TD[2] = &Cp;
  varnames_TD[3] = "ViscosityDyn";
  val_vars_TD[3] = &Mu;
  varnames_TD[4] = "Conductivity";
  val_vars_TD[4] = &Kt;
  varnames_TD[5] = "DiffusionCoefficient";
  val_vars_TD[5] = &mass_diffusivity;
  varnames_TD[6] = "MolarWeightMix";
  val_vars_TD[6] = &molar_weight;

  /*--- Source term variables ---*/
  varnames_Sources.resize(n_table_sources);
  val_vars_Sources.resize(n_table_sources);
  
  for(size_t iSource=0;iSource<n_table_sources; iSource++){
    varnames_Sources[iSource] = table_source_names[iSource];
    val_vars_Sources[iSource] = &table_sources[iSource];
  }

  /*--- Passive lookups ---*/
  varnames_LookUp.resize(n_lookups);
  val_vars_LookUp.resize(n_lookups);

  for(size_t iLookup=0; iLookup<n_lookups; iLookup++){
    varnames_LookUp[iLookup] = table_lookup_names[iLookup];
    val_vars_LookUp[iLookup] = &lookup_scalar[iLookup];
  }
}