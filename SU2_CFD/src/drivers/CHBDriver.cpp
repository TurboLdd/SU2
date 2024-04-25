/*!
 * \file CDriver.cpp
 * \brief The main subroutines for driving single or multi-zone problems.
 * \author T. Economon, H. Kline, R. Sanchez, F. Palacios
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/drivers/CHBDriver.hpp"
#include "../../include/definition_structure.hpp"
#include "../../include/output/COutput.hpp"
#include "../../include/iteration/CIteration.hpp"
CHBDriver::CHBDriver(char* confFile,
    unsigned short val_nZone,
    SU2_Comm MPICommunicator) : CFluidDriver(confFile,
        val_nZone,
        MPICommunicator) {
  unsigned short kInst;

  nInstHB = nInst[ZONE_0];

  /*--- allocate dynamic memory for the Harmonic Balance operator ---*/
  D = new su2double*[nInstHB];
  Sr = new su2double*[nInstHB];
  Sl = new su2double*[nInstHB];
  for (kInst = 0; kInst < nInstHB; kInst++) {
    D[kInst] = new su2double[nInstHB];
    Sr[kInst] = new su2double[nInstHB];
    Sl[kInst] = new su2double[nInstHB];
  }
}

CHBDriver::~CHBDriver() {

  unsigned short kInst;

  /*--- delete dynamic memory for the Harmonic Balance operator ---*/
  for (kInst = 0; kInst < nInstHB; kInst++) {
    delete[] D[kInst];
    delete[] Sr[kInst];
    delete[] Sl[kInst];
  }
  delete [] D,Sl,Sr;
}


void CHBDriver::Run() {

  /*--- Run a single iteration of a Harmonic Balance problem. Preprocess all
   all zones before beginning the iteration. ---*/

  for (iInst = 0; iInst < nInstHB; iInst++){
    iteration_container[ZONE_0][iInst]->Preprocess(output_container[ZONE_0], integration_container, geometry_container,
        solver_container, numerics_container, config_container,
        surface_movement, grid_movement, FFDBox, ZONE_0, iInst);}

  for (iInst = 0; iInst < nInstHB; iInst++)
    iteration_container[ZONE_0][iInst]->Iterate(output_container[ZONE_0], integration_container, geometry_container,
        solver_container, numerics_container, config_container,
        surface_movement, grid_movement, FFDBox, ZONE_0, iInst);

  for (iInst = 0; iInst < nInstHB; iInst++)
    iteration_container[ZONE_0][iInst]->Monitor(output_container[ZONE_0], integration_container, geometry_container,
        solver_container, numerics_container, config_container,
        surface_movement, grid_movement, FFDBox, ZONE_0, iInst);

}

void CHBDriver::Update() {


  //1, imarker
  //2, iInst
  //3,ivertex
  //4, iVar
  for (iInst = 0; iInst < nInstHB; iInst++) {
    /*--- Compute the harmonic balance terms across all zones ---*/
    SetHarmonicBalance(iInst);

  }

  /*--- Precondition the harmonic balance source terms ---*/
  if (config_container[ZONE_0]->GetHB_Precondition() == YES) {
    StabilizeHarmonicBalance();

  }
  if (config_container[ZONE_0]->GetHB_PhaseLag() == YES) {
    
      SetHarmonicBalancePhaseLag();

    
  }

  for (iInst = 0; iInst < nInstHB; iInst++) {

    /*--- Update the harmonic balance terms across all zones ---*/
    iteration_container[ZONE_0][iInst]->Update(output_container[ZONE_0], integration_container, geometry_container,
        solver_container, numerics_container, config_container,
        surface_movement, grid_movement, FFDBox, ZONE_0, iInst);

  }

}

void CHBDriver::SetHarmonicBalance(unsigned short iInst) {

  unsigned short iVar, jInst, iMGlevel;
  unsigned short nVar = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetnVar();
  unsigned long iPoint;
  bool implicit = (config_container[ZONE_0]->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool adjoint = (config_container[ZONE_0]->GetContinuous_Adjoint());
  if (adjoint) {
    implicit = (config_container[ZONE_0]->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  }

  unsigned long InnerIter = config_container[ZONE_0]->GetInnerIter();

  /*--- Retrieve values from the config file ---*/
  auto *U = new su2double[nVar];
  auto *U_old = new su2double[nVar];
  auto *Psi = new su2double[nVar];
  auto *Psi_old = new su2double[nVar];
  auto *Source = new su2double[nVar];
  su2double deltaU, deltaPsi;

  /*--- Compute period of oscillation ---*/
  su2double period = config_container[ZONE_0]->GetHarmonicBalance_Period();

  /*--- Non-dimensionalize the input period, if necessary.  */
  period /= config_container[ZONE_0]->GetTime_Ref();

  if (InnerIter == 0)
    ComputeHBOperator();

  /*--- Compute various source terms for explicit direct, implicit direct, and adjoint problems ---*/
  /*--- Loop over all grid levels ---*/
  for (iMGlevel = 0; iMGlevel <= config_container[ZONE_0]->GetnMGLevels(); iMGlevel++) {

    /*--- Loop over each node in the volume mesh ---*/
    SU2_OMP_FOR_DYN(256)//TODO: Check if this is the best number of threads
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][iInst][iMGlevel]->GetnPoint(); iPoint++) {

      for (iVar = 0; iVar < nVar; iVar++) {
        Source[iVar] = 0.0;
      }

      /*--- Step across the columns ---*/
      for (jInst = 0; jInst < nInstHB; jInst++) {

        /*--- Retrieve solution at this node in current zone ---*/
        for (iVar = 0; iVar < nVar; iVar++) {

          if (!adjoint) {
            U[iVar] = solver_container[ZONE_0][jInst][iMGlevel][FLOW_SOL]->GetNodes()->GetSolution(iPoint, iVar);
            Source[iVar] += U[iVar]*D[iInst][jInst];

            if (implicit) {
              U_old[iVar] = solver_container[ZONE_0][jInst][iMGlevel][FLOW_SOL]->GetNodes()->GetSolution_Old(iPoint, iVar);
              deltaU = U[iVar] - U_old[iVar];
              Source[iVar] += deltaU*D[iInst][jInst];
            }

          }

          else {
            Psi[iVar] = solver_container[ZONE_0][jInst][iMGlevel][ADJFLOW_SOL]->GetNodes()->GetSolution(iPoint, iVar);
            Source[iVar] += Psi[iVar]*D[jInst][iInst];

            if (implicit) {
              Psi_old[iVar] = solver_container[ZONE_0][jInst][iMGlevel][ADJFLOW_SOL]->GetNodes()->GetSolution_Old(iPoint, iVar);
              deltaPsi = Psi[iVar] - Psi_old[iVar];
              Source[iVar] += deltaPsi*D[jInst][iInst];
            }
          }
        }

        /*--- Store sources for current row ---*/
        for (iVar = 0; iVar < nVar; iVar++) {
          if (!adjoint) {
            solver_container[ZONE_0][iInst][iMGlevel][FLOW_SOL]->GetNodes()->SetHarmonicBalance_Source(iPoint, iVar, Source[iVar]);
          }
          else {
            solver_container[ZONE_0][iInst][iMGlevel][ADJFLOW_SOL]->GetNodes()->SetHarmonicBalance_Source(iPoint, iVar, Source[iVar]);
          }
        }

      }
    }
    END_SU2_OMP_FOR


  }

  /*--- Source term for a turbulence model ---*/
  if (config_container[ZONE_0]->GetKind_Solver() == MAIN_SOLVER::RANS) {

    /*--- Extra variables needed if we have a turbulence model. ---*/
    unsigned short nVar_Turb = solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->GetnVar();
    auto *U_Turb = new su2double[nVar_Turb];
    auto *Source_Turb = new su2double[nVar_Turb];

    /*--- Loop over only the finest mesh level (turbulence is always solved
     on the original grid only). ---*/
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint(); iPoint++) {
      for (iVar = 0; iVar < nVar_Turb; iVar++) Source_Turb[iVar] = 0.0;
      for (jInst = 0; jInst < nInstHB; jInst++) {

        /*--- Retrieve solution at this node in current zone ---*/
        for (iVar = 0; iVar < nVar_Turb; iVar++) {
          U_Turb[iVar] = solver_container[ZONE_0][jInst][MESH_0][TURB_SOL]->GetNodes()->GetSolution(iPoint, iVar);
          Source_Turb[iVar] += U_Turb[iVar]*D[iInst][jInst];
        }
      }

      /*--- Store sources for current iZone ---*/
      for (iVar = 0; iVar < nVar_Turb; iVar++)
        solver_container[ZONE_0][iInst][MESH_0][TURB_SOL]->GetNodes()->SetHarmonicBalance_Source(iPoint, iVar, Source_Turb[iVar]);
    }

    delete [] U_Turb;
    delete [] Source_Turb;
  }

  delete [] Source;
  delete [] U;
  delete [] U_old;
  delete [] Psi;
  delete [] Psi_old;

}
/*currently adjoint phaselag is not supported*/
void CHBDriver::SetHarmonicBalancePhaseLag() {
  unsigned short iInst, jInst;
  unsigned short iVar, jVar, iMGlevel = 0;
  unsigned short nVar = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetnVar();
  unsigned long iVertex, iPoint, iPeriodic;
  unsigned short iMarker;
  bool implicit = (config_container[ZONE_0]->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool adjoint = (config_container[ZONE_0]->GetContinuous_Adjoint());
  bool ransturb = (config_container[ZONE_0]->GetKind_Solver() == MAIN_SOLVER::RANS);
  if (adjoint) {
    SU2_MPI::Error("Turn off the adjoint solve when you want to use Phase Lag BC.", CURRENT_FUNCTION);
    // implicit = (config_container[ZONE_0]->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  }
  unsigned short nVar_Turb=1;//very ugly!!!here i want to loop the vertex only once to deduce the memory reading,so i defined the pointer here
  if (ransturb) {
    nVar_Turb = solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->GetnVar();
  }

  /*--- Retrieve values from the config file ---*/
  auto* U = new su2double[nVar];
  auto** dU = new su2double*[nVar];
  auto* U2U = new su2double[nVar];      // right to left
  auto** dU2dU = new su2double*[nVar];  // implicit right to left

  for (iVar = 0; iVar < nVar; iVar++) {
    dU2dU[iVar] = new su2double[nVar];
    dU[iVar] = new su2double[nVar];
  }
  auto* U_Turb = new su2double[nVar_Turb];
  auto* U2U_Turb = new su2double[nVar_Turb];
  auto** dU_Turb = new su2double*[nVar_Turb];
  auto** dU2dU_Turb = new su2double*[nVar_Turb];
  for (iVar = 0; iVar < nVar_Turb; iVar++) {
    dU_Turb[iVar] = new su2double[nVar_Turb];
    dU2dU_Turb[iVar] = new su2double[nVar_Turb];
  }
  // auto* Psi = new su2double[nVar];// adjoint variable
  // auto* Psi_old = new su2double[nVar];//implicit adjoint variable
  //  initialize the array

  // su2double deltaU, deltaPsi;




  //

  /*--- Compute various bc  for explicit direct, implicit direct, and adjoint problems ---*/
  /*--- Loop over all grid levels ---*/
  // periodic BC should be solved on finest mesh
  // for (iMGlevel = 0; iMGlevel <= config_container[ZONE_0]->GetnMGLevels(); iMGlevel++) {
  /*--- Step1: determined the periodic pair to be matched ---*/
  unsigned short nPeriodic = config_container[ZONE_0]->GetnMarker_Periodic() / 2;
  // 1,find the periodic pair
  for (unsigned short val_periodic = 1; val_periodic <= nPeriodic; val_periodic++) {
    for (iMarker = 0; iMarker < config_container[ZONE_0]->GetnMarker_All(); iMarker++) {
      if (config_container[ZONE_0]->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
        iPeriodic = config_container[ZONE_0]->GetMarker_All_PerBound(iMarker);
        if ((iPeriodic == val_periodic) || (iPeriodic == val_periodic + nPeriodic)) {
          /*--- Loop over each node in the periodic mesh ---*/
          // 2,loop the node
          for (iInst = 0; iInst < nInstHB; iInst++) {
            SU2_OMP_FOR_DYN(256)  // TODO: Check if this is the best number of threads
            for (iVertex = 0; iVertex < geometry_container[ZONE_0][iInst][iMGlevel]->nVertex[iMarker]; iVertex++) {
              iPoint = geometry_container[ZONE_0][iInst][iMGlevel]->vertex[iMarker][iVertex]->GetNode();
              /*--- Step across the columns ---*/
              // initialize the temp array for each iPoint
              for (iVar = 0; iVar < nVar; iVar++) {
                U2U[iVar] = 0.0;

                if (implicit) {
                  for (jVar = 0; jVar < nVar; jVar++) {
                    dU2dU[iVar][jVar] = 0.0;
                  }
                }
              }
              if (ransturb) {
                for (iVar = 0; iVar < nVar_Turb; iVar++) {
                  U_Turb[iVar] = 0.0;
                  if (implicit) {
                    for (jVar = 0; jVar < nVar_Turb; jVar++) {
                      dU2dU_Turb[iVar][jVar] = 0.0;
                    }
                  }
                }
              }

              // 3,hb cross
              for (jInst = 0; jInst < nInstHB; jInst++) {
                /*--- Retrieve solution at this node in current zone ---*/
                // 4,loop variable
                for (iVar = 0; iVar < nVar; iVar++) {
                  if (!adjoint) {
                    U[iVar] = solver_container[ZONE_0][jInst][iMGlevel][FLOW_SOL]->LinSysRes(iPoint, iVar);
                    if (iPeriodic == val_periodic) {
                      U2U[iVar] += Sl[iInst][jInst] * U[iVar];
                    } else if (iPeriodic == val_periodic + nPeriodic) {
                      U2U[iVar] += Sr[iInst][jInst] * U[iVar];
                    }
                    // get jacobian
                    if (implicit) {
                      for (jVar = 0; jVar < nVar; jVar++) {
                        dU[iVar][jVar] = solver_container[ZONE_0][jInst][iMGlevel][FLOW_SOL]->Jacobian.GetBlock(
                            iPoint, iPoint, iVar, jVar);
                        if (iPeriodic == val_periodic) {
                          dU2dU[iVar][jVar] += Sl[iInst][jInst] * dU[iVar][jVar];
                        } else if (iPeriodic == val_periodic + nPeriodic) {
                          dU2dU[iVar][jVar] += Sr[iInst][jInst] * dU[iVar][jVar];
                        }
                      }
                    }
                  }
                  /*--- Phase Lag for a turbulence model ---*/
                  if (ransturb) {
                    /*--- Extra variables needed if we have a turbulence model. ---*/

                    /*--- Retrieve solution at this node in current zone ---*/
                    for (iVar = 0; iVar < nVar_Turb; iVar++) {
                      U_Turb[iVar] =
                          solver_container[ZONE_0][jInst][MESH_0][TURB_SOL]->GetNodes()->GetSolution(iPoint, iVar);
                      if (iPeriodic == val_periodic) {
                        U2U_Turb[iVar] += U_Turb[iVar] * Sl[iInst][jInst];
                      } else if (iPeriodic == val_periodic + nPeriodic) {
                        U2U_Turb[iVar] += U_Turb[iVar] * Sr[iInst][jInst];
                      }
                      if (implicit) {
                        for (jVar = 0; jVar < nVar_Turb; jVar++) {
                          dU_Turb[iVar][jVar] = solver_container[ZONE_0][jInst][iMGlevel][TURB_SOL]->Jacobian.GetBlock(
                              iPoint, iPoint, iVar, jVar);
                          if (iPeriodic == val_periodic) {
                            dU2dU_Turb[iVar][jVar] += Sl[iInst][jInst] * dU_Turb[iVar][jVar];
                          } else if (iPeriodic == val_periodic + nPeriodic) {
                            dU2dU_Turb[iVar][jVar] += Sr[iInst][jInst] * dU_Turb[iVar][jVar];
                          }
                        }
                      }
                    }
                  }

                  // else {
                  //  Psi[iVar] = solver_container[ZONE_0][jInst][iMGlevel][ADJFLOW_SOL]->LinSysRes(iPoint, iVar);

                  /*if (implicit) {
                    Psi_old[iVar] =
                  solver_container[ZONE_0][jInst][iMGlevel][ADJFLOW_SOL]->GetNodes()->GetSolution_Old( iPoint, iVar);
                    deltaPsi = Psi[iVar] - Psi_old[iVar];
                    // need to add the matrix production
                  }*/
                  //}
                }

                /*--- Store sources for current row ---*/
                for (iVar = 0; iVar < nVar; iVar++) {
                  if (!adjoint) {
                    if (iPeriodic == val_periodic) {
                      solver_container[ZONE_0][iInst][iMGlevel][FLOW_SOL]->LinSysRes.AddBlock(iPoint, U2U);
                      if (implicit) {
                        solver_container[ZONE_0][iInst][iMGlevel][FLOW_SOL]->Jacobian.AddBlock2Diag(iPoint, dU2dU);
                      }
                    } else if (iPeriodic == val_periodic + nPeriodic) {
                      solver_container[ZONE_0][iInst][iMGlevel][FLOW_SOL]->LinSysRes.AddBlock(iPoint, U2U);
                      if (implicit) {
                        solver_container[ZONE_0][iInst][iMGlevel][FLOW_SOL]->Jacobian.AddBlock2Diag(iPoint, dU2dU);
                      }
                    }
                  }
                  if (ransturb) {
                    if (iPeriodic == val_periodic) {
                      solver_container[ZONE_0][iInst][iMGlevel][TURB_SOL]->LinSysRes.SubtractBlock(iPoint, U2U);
                      if (implicit) {
                        solver_container[ZONE_0][iInst][iMGlevel][TURB_SOL]->Jacobian.AddBlock2Diag(iPoint, dU2dU);
                      }
                    } else if (iPeriodic == val_periodic + nPeriodic) {
                      solver_container[ZONE_0][iInst][iMGlevel][TURB_SOL]->LinSysRes.SubtractBlock(iPoint, U2U_Turb);
                      if (implicit) {
                        solver_container[ZONE_0][iInst][iMGlevel][TURB_SOL]->Jacobian.AddBlock2Diag(iPoint, dU2dU_Turb);
                      }
                    }
                  }
                  /*else {
                    if (iPeriodic == val_periodic) {
                      solver_container[ZONE_0][iInst][iMGlevel][ADJFLOW_SOL]->SetCharacPrimVarPhaseLagR(
                          iMarker, iVertex, iVar, Ur2Ul[iVar]);
                    } else if (iPeriodic == val_periodic + nPeriodic) {
                      solver_container[ZONE_0][iInst][iMGlevel][ADJFLOW_SOL]->SetCharacPrimVarPhaseLagL(
                          iMarker, iVertex, iVar, Ul2Ur[iVar]);
                    }
                  }*/
                }
              }
            }
            END_SU2_OMP_FOR
          }
        }
      }
    }
  }
  //}

  for (iVar = 0; iVar < nVar; iVar++) {
    delete[] dU2dU[iVar];
    delete[] dU[iVar];
  }
  for (iVar = 0; iVar < nVar_Turb; iVar++) {
    delete[] dU2dU_Turb[iVar];
    delete[] dU_Turb[iVar];
  }
  delete[] U, dU, U2U, dU2dU;
  delete[] U_Turb,U2U_Turb,dU_Turb,dU2dU_Turb;


  // delete[] Psi;
  // delete[] Psi_old;
}

void CHBDriver::StabilizeHarmonicBalance() {

  unsigned short i, j, k, iVar, iInst, jInst, iMGlevel;
  unsigned short nVar = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetnVar();
  unsigned long iPoint;
  bool adjoint = (config_container[ZONE_0]->GetContinuous_Adjoint());

  /*--- Retrieve values from the config file ---*/
  auto *Source     = new su2double[nInstHB];
  auto *Source_old = new su2double[nInstHB];
  su2double Delta;

  auto **Pinv     = new su2double*[nInstHB];
  auto **P        = new su2double*[nInstHB];
  for (iInst = 0; iInst < nInstHB; iInst++) {
    Pinv[iInst]       = new su2double[nInstHB];
    P[iInst]          = new su2double[nInstHB];
  }

  /*--- Loop over all grid levels ---*/
  for (iMGlevel = 0; iMGlevel <= config_container[ZONE_0]->GetnMGLevels(); iMGlevel++) {

    /*--- Loop over each node in the volume mesh ---*/
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][INST_0][iMGlevel]->GetnPoint(); iPoint++) {

      /*--- Get time step for current node ---*/
      Delta = solver_container[ZONE_0][INST_0][iMGlevel][FLOW_SOL]->GetNodes()->GetDelta_Time(iPoint);

      /*--- Setup stabilization matrix for this node ---*/
      for (iInst = 0; iInst < nInstHB; iInst++) {
        for (jInst = 0; jInst < nInstHB; jInst++) {
          if (jInst == iInst ) {
            Pinv[iInst][jInst] = 1.0 + Delta*D[iInst][jInst];
          }
          else {
            Pinv[iInst][jInst] = Delta*D[iInst][jInst];
          }
        }
      }

      /*--- Invert stabilization matrix Pinv with Gauss elimination---*/

      /*--  A temporary matrix to hold the inverse, dynamically allocated ---*/
      auto **temp = new su2double*[nInstHB];
      for (i = 0; i < nInstHB; i++) {
        temp[i] = new su2double[2 * nInstHB];
      }

      /*---  Copy the desired matrix into the temporary matrix ---*/
      for (i = 0; i < nInstHB; i++) {
        for (j = 0; j < nInstHB; j++) {
          temp[i][j] = Pinv[i][j];
          temp[i][nInstHB + j] = 0;
        }
        temp[i][nInstHB + i] = 1;
      }

      su2double max_val;
      unsigned short max_idx;

      /*---  Pivot each column such that the largest number possible divides the other rows  ---*/
      for (k = 0; k < nInstHB - 1; k++) {
        max_idx = k;
        max_val = abs(temp[k][k]);
        /*---  Find the largest value (pivot) in the column  ---*/
        for (j = k; j < nInstHB; j++) {
          if (abs(temp[j][k]) > max_val) {
            max_idx = j;
            max_val = abs(temp[j][k]);
          }
        }

        /*---  Move the row with the highest value up  ---*/
        for (j = 0; j < (nInstHB * 2); j++) {
          su2double d = temp[k][j];
          temp[k][j] = temp[max_idx][j];
          temp[max_idx][j] = d;
        }
        /*---  Subtract the moved row from all other rows ---*/
        for (i = k + 1; i < nInstHB; i++) {
          su2double c = temp[i][k] / temp[k][k];
          for (j = 0; j < (nInstHB * 2); j++) {
            temp[i][j] = temp[i][j] - temp[k][j] * c;
          }
        }
      }

      /*---  Back-substitution  ---*/
      for (k = nInstHB - 1; k > 0; k--) {
        if (temp[k][k] != su2double(0.0)) {
          for (int i = k - 1; i > -1; i--) {
            su2double c = temp[i][k] / temp[k][k];
            for (j = 0; j < (nInstHB * 2); j++) {
              temp[i][j] = temp[i][j] - temp[k][j] * c;
            }
          }
        }
      }

      /*---  Normalize the inverse  ---*/
      for (i = 0; i < nInstHB; i++) {
        su2double c = temp[i][i];
        for (j = 0; j < nInstHB; j++) {
          temp[i][j + nInstHB] = temp[i][j + nInstHB] / c;
        }
      }

      /*---  Copy the inverse back to the main program flow ---*/
      for (i = 0; i < nInstHB; i++) {
        for (j = 0; j < nInstHB; j++) {
          P[i][j] = temp[i][j + nInstHB];
        }
      }

      /*---  Delete dynamic template  ---*/
      for (iInst = 0; iInst < nInstHB; iInst++) {
        delete[] temp[iInst];
      }
      delete[] temp;

      /*--- Loop through variables to precondition ---*/
      for (iVar = 0; iVar < nVar; iVar++) {

        /*--- Get current source terms (not yet preconditioned) and zero source array to prepare preconditioning ---*/
        for (iInst = 0; iInst < nInstHB; iInst++) {
          Source_old[iInst] = solver_container[ZONE_0][iInst][iMGlevel][FLOW_SOL]->GetNodes()->GetHarmonicBalance_Source(iPoint, iVar);
          Source[iInst] = 0;
        }

        /*--- Step through columns ---*/
        for (iInst = 0; iInst < nInstHB; iInst++) {
          for (jInst = 0; jInst < nInstHB; jInst++) {
            Source[iInst] += P[iInst][jInst]*Source_old[jInst];
          }

          /*--- Store updated source terms for current node ---*/
          if (!adjoint) {
            solver_container[ZONE_0][iInst][iMGlevel][FLOW_SOL]->GetNodes()->SetHarmonicBalance_Source(iPoint, iVar, Source[iInst]);
          }
          else {
            solver_container[ZONE_0][iInst][iMGlevel][ADJFLOW_SOL]->GetNodes()->SetHarmonicBalance_Source(iPoint, iVar, Source[iInst]);
          }
        }

      }
    }
  }

  /*--- Deallocate dynamic memory ---*/
  for (iInst = 0; iInst < nInstHB; iInst++){
    delete [] P[iInst];
    delete [] Pinv[iInst];
  }
  delete [] P;
  delete [] Pinv;
  delete [] Source;
  delete [] Source_old;

}

void CHBDriver::ComputeHBOperator() {

  const   complex<su2double> J(0.0,1.0);
  unsigned short i, j, k, iInst;

  auto *Omega_HB       = new su2double[nInstHB];
  auto *TimeLag_HB       = new su2double[nInstHB];
  auto **E    = new complex<su2double>*[nInstHB];
  auto **Einv = new complex<su2double>*[nInstHB];
  auto **DD   = new complex<su2double>*[nInstHB];
  auto **SSR   = new complex<su2double>*[nInstHB];
  auto **SSL   = new complex<su2double>*[nInstHB];
  for (iInst = 0; iInst < nInstHB; iInst++) {
    E[iInst]    = new complex<su2double>[nInstHB];
    Einv[iInst] = new complex<su2double>[nInstHB];
    DD[iInst]   = new complex<su2double>[nInstHB];
    SSR[iInst]   = new complex<su2double>[nInstHB];
    SSL[iInst]   = new complex<su2double>[nInstHB];
  }

  /*--- Get simualation period from config file ---*/
  su2double Period = config_container[ZONE_0]->GetHarmonicBalance_Period();

  /*--- Non-dimensionalize the input period, if necessary.      */
  Period /= config_container[ZONE_0]->GetTime_Ref();

  /*--- Build the array containing the selected frequencies to solve ---*/
  for (iInst = 0; iInst < nInstHB; iInst++) {
    Omega_HB[iInst]  = config_container[ZONE_0]->GetOmega_HB()[iInst];
    Omega_HB[iInst] /= config_container[ZONE_0]->GetOmega_Ref(); //TODO: check
    TimeLag_HB[iInst]  = config_container[ZONE_0]->GetPhaseLag_HB()[iInst];// Now, only uniform time sampling
  }

  /*--- Build the diagonal matrix of the frequencies DD and Phase Shifted SS---*/
  for (i = 0; i < nInstHB; i++) {
    for (k = 0; k < nInstHB; k++) {
      if (k == i ) {
        DD[i][k] = J*Omega_HB[k];
        if(config_container[ZONE_0]->GetHB_PhaseLag() == YES){
          SSR[i][k] = complex<su2double>(cos(TimeLag_HB[k])) + J*complex<su2double>(sin(TimeLag_HB[k]));
          SSL[i][k] = complex<su2double>(cos(TimeLag_HB[k])) - J*complex<su2double>(sin(TimeLag_HB[k]));
        }
       
      }
    }
  }


  /*--- Build the harmonic balance inverse matrix ---*/
  for (i = 0; i < nInstHB; i++) {
    for (k = 0; k < nInstHB; k++) {
      Einv[i][k] = complex<su2double>(cos(Omega_HB[k]*(i*Period/nInstHB))) + J*complex<su2double>(sin(Omega_HB[k]*(i*Period/nInstHB)));
    }
  }

  /*---  Invert inverse harmonic balance Einv with Gauss elimination ---*/

  /*--  A temporary matrix to hold the inverse, dynamically allocated ---*/
  auto **temp = new complex<su2double>*[nInstHB];
  for (i = 0; i < nInstHB; i++) {
    temp[i] = new complex<su2double>[2 * nInstHB];
  }

  /*---  Copy the desired matrix into the temporary matrix ---*/
  for (i = 0; i < nInstHB; i++) {
    for (j = 0; j < nInstHB; j++) {
      temp[i][j] = Einv[i][j];
      temp[i][nInstHB + j] = 0;
    }
    temp[i][nInstHB + i] = 1;
  }

  su2double max_val;
  unsigned short max_idx;

  /*---  Pivot each column such that the largest number possible divides the other rows  ---*/
  for (k = 0; k < nInstHB - 1; k++) {
    max_idx = k;
    max_val = abs(temp[k][k]);
    /*---  Find the largest value (pivot) in the column  ---*/
    for (j = k; j < nInstHB; j++) {
      if (abs(temp[j][k]) > max_val) {
        max_idx = j;
        max_val = abs(temp[j][k]);
      }
    }
    /*---  Move the row with the highest value up  ---*/
    for (j = 0; j < (nInstHB * 2); j++) {
      complex<su2double> d = temp[k][j];
      temp[k][j] = temp[max_idx][j];
      temp[max_idx][j] = d;
    }
    /*---  Subtract the moved row from all other rows ---*/
    for (i = k + 1; i < nInstHB; i++) {
      complex<su2double> c = temp[i][k] / temp[k][k];
      for (j = 0; j < (nInstHB * 2); j++) {
        temp[i][j] = temp[i][j] - temp[k][j] * c;
      }
    }
  }
  /*---  Back-substitution  ---*/
  for (k = nInstHB - 1; k > 0; k--) {
    if (temp[k][k] != complex<su2double>(0.0)) {
      for (int i = k - 1; i > -1; i--) {
        complex<su2double> c = temp[i][k] / temp[k][k];
        for (j = 0; j < (nInstHB * 2); j++) {
          temp[i][j] = temp[i][j] - temp[k][j] * c;
        }
      }
    }
  }
  /*---  Normalize the inverse  ---*/
  for (i = 0; i < nInstHB; i++) {
    complex<su2double> c = temp[i][i];
    for (j = 0; j < nInstHB; j++) {
      temp[i][j + nInstHB] = temp[i][j + nInstHB] / c;
    }
  }
  /*---  Copy the inverse back to the main program flow ---*/
  for (i = 0; i < nInstHB; i++) {
    for (j = 0; j < nInstHB; j++) {
      E[i][j] = temp[i][j + nInstHB];
    }
  }
  /*---  Delete dynamic template  ---*/
  for (i = 0; i < nInstHB; i++) {
    delete[] temp[i];
  }
  delete[] temp;


  /*---  Temporary matrix for performing product  ---*/
  auto **Temp    = new complex<su2double>*[nInstHB];
  auto **SRTemp    = new complex<su2double>*[nInstHB];
  auto **SLTemp    = new complex<su2double>*[nInstHB];

  /*---  Temporary complex HB operator  ---*/
  auto **Dcpx    = new complex<su2double>*[nInstHB];
  auto **SRcpx    = new complex<su2double>*[nInstHB];
  auto **SLcpx    = new complex<su2double>*[nInstHB];

  for (iInst = 0; iInst < nInstHB; iInst++){
    Temp[iInst]    = new complex<su2double>[nInstHB];
    SRTemp[iInst]    = new complex<su2double>[nInstHB];
    SLTemp[iInst]    = new complex<su2double>[nInstHB];
    Dcpx[iInst]   = new complex<su2double>[nInstHB];
    SRcpx[iInst]   = new complex<su2double>[nInstHB];
    SLcpx[iInst]   = new complex<su2double>[nInstHB];
  }


  /*---  Calculation of the HB operator matrix ---*/
  for (int row = 0; row < nInstHB; row++) {
    for (int col = 0; col < nInstHB; col++) {
      for (int inner = 0; inner < nInstHB; inner++) {
        Temp[row][col] += Einv[row][inner] * DD[inner][col];
        if(config_container[ZONE_0]->GetHB_PhaseLag() == YES){
        SRTemp[row][col] += Einv[row][inner] * SSR[inner][col];
        SLTemp[row][col] += Einv[row][inner] * SSL[inner][col];
        }
      }
    }
  }

  unsigned short row, col, inner;

  for (row = 0; row < nInstHB; row++) {
    for (col = 0; col < nInstHB; col++) {
      for (inner = 0; inner < nInstHB; inner++) {
        Dcpx[row][col] += Temp[row][inner] * E[inner][col];
        if(config_container[ZONE_0]->GetHB_PhaseLag() == YES){
        SRcpx[row][col] += SRTemp[row][inner] * E[inner][col];
        SLcpx[row][col] += SRTemp[row][inner] * E[inner][col];
        }
      }
    }
  }

  /*---  Take just the real part of the HB operator matrix ---*/
  for (i = 0; i < nInstHB; i++) {
    for (k = 0; k < nInstHB; k++) {
      D[i][k] = real(Dcpx[i][k]);
      if(config_container[ZONE_0]->GetHB_PhaseLag() == YES){
      Sr[i][k] = real(SRcpx[i][k]);
      Sl[i][k] = real(SLcpx[i][k]);
      }
    }
  }

  /*--- Deallocate dynamic memory ---*/
  for (iInst = 0; iInst < nInstHB; iInst++){
    delete [] E[iInst];
    delete [] Einv[iInst];
    delete [] DD[iInst];
    delete [] SSR[iInst];
    delete [] SSL[iInst];
    delete [] Temp[iInst];
    delete [] SRTemp[iInst];
    delete [] SLTemp[iInst];
    delete [] Dcpx[iInst];
    delete [] SRcpx[iInst];
    delete [] SLcpx[iInst];
  }
  delete [] E;
  delete [] Einv;
  delete [] DD;
  delete [] SSR;
  delete [] SSL;
  delete [] Temp;
  delete [] SRTemp;
  delete [] SLTemp;
  delete [] Dcpx;
  delete [] SRcpx;
  delete [] SLcpx;
  delete [] Omega_HB,TimeLag_HB;

}
void CHBDriver::Postprocess(unsigned short iInst) {

  iteration_container[ZONE_0][iInst]->Postprocess(output_container[ZONE_0], integration_container, geometry_container, solver_container,
      numerics_container, config_container, surface_movement, grid_movement, FFDBox, ZONE_0, iInst);

  /*--- A corrector step can help preventing numerical instabilities ---*/

  if (config_container[ZONE_0]->GetRelaxation())
    iteration_container[ZONE_0][iInst]->Relaxation(output_container[ZONE_0], integration_container, geometry_container, solver_container,
        numerics_container, config_container, surface_movement, grid_movement, FFDBox, ZONE_0, iInst);

}