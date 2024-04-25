/*!
 * \class CHBDriver
 * \ingroup Drivers
 * \brief Class for driving an iteration of Harmonic Balance (HB) method problem using multiple time zones.
 * \author T. Economon
 */
#pragma once


#include "CDriverBase.hpp"
#include "CDriver.hpp"


#include <cassert>

#ifdef VTUNEPROF
#include <ittnotify.h>
#endif
#include <cfenv>
using namespace std;
class CHBDriver : public CFluidDriver {
 private:
  unsigned short nInstHB;
  su2double** D; /*!< \brief Harmonic Balance operator. */
  su2double** Sr; /*!< \brief Harmonic Balance operator for phase lag bc. */
  su2double** Sl; /*!< \brief Harmonic Balance operator for phase lag bc. */

  /*!
   * \brief Computation and storage of the Harmonic Balance method source terms.
   * \author T. Economon, K. Naik
   * \param[in] iZone - Current zone number.
   */
  void SetHarmonicBalance(unsigned short iZone);

    /*!
   * \brief Computation and storage of the Harmonic Balance method phase lag bc.
   * \author Xuedong Zheng
   * \param[in] iZone - Current zone number.
   */
  void SetHarmonicBalancePhaseLag();
  /*!
   * \brief Precondition Harmonic Balance source term for stability
   * \author J. Howison
   */
  void StabilizeHarmonicBalance();

  /*!
   * \brief Computation of the Harmonic Balance operator matrix for harmonic balance.
   * \author A. Rubino, S. Nimmagadda
   */
  void ComputeHBOperator();

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] confFile - Configuration file name.
   * \param[in] val_nZone - Total number of zones.
   * \param[in] MPICommunicator - MPI communicator for SU2.
   */
  CHBDriver(char* confFile, unsigned short val_nZone, SU2_Comm MPICommunicator);

  /*!
   * \brief Destructor of the class.
   */
  ~CHBDriver(void) override;

  /*!
   * \brief Run a single iteration of a Harmonic Balance problem.
   */
  void Run() override;

  /*!
   * \brief Update the solution for the Harmonic Balance.
   */
  void Update() override;

  inline unsigned short GetnInstHB() const { return nInstHB; }

  /*!
   * \brief Set the mesh displacements of a marker vertex.
   * \note This can be the input of the flow solver in an FSI setting.
   * \param[in] iMarker - Marker index.
   * \param[in] iVertex - Marker vertex index.
   * \param[in] values - Node displacements (nDim).
   */
  inline void SetMarkerCustomDisplacement(unsigned short iMarker, unsigned long iVertex, vector<passivedouble> values,
                                          unsigned short iInst) {
    const auto iPoint = GetMarkerNode(iMarker, iVertex, iInst);          //
    auto* nodes = GetSolverAndCheckMarker(MESH_SOL, iInst)->GetNodes();  //

    for (auto iDim = 0u; iDim < GetNumberDimensions(); iDim++) {  //
      nodes->SetBound_Disp(iPoint, iDim, values[iDim]);           //
    }
  }

  void Postprocess() ;
  
  unsigned long GetMarkerNode(unsigned short iMarker, unsigned long iVertex, unsigned short iInst) const {
    if (iVertex >= GetNumberMarkerNodes(iMarker, iInst)) {  //
      SU2_MPI::Error("Vertex index exceeds marker size.", CURRENT_FUNCTION);
    }
    return geometry_container[ZONE_0][iInst][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  }
  unsigned long GetNumberMarkerNodes(unsigned short iMarker, unsigned short iInst) const {  //
    if (iMarker >= GetNumberMarkers(ZONE_0)) {                                              //
      SU2_MPI::Error("Marker index exceeds size.", CURRENT_FUNCTION);
    }
    return geometry_container[ZONE_0][iInst][MESH_0]->GetnVertex(iMarker);
  }
  unsigned short GetNumberMarkers(unsigned short iZone) const { return config_container[iZone]->GetnMarker_All(); }

  inline CSolver* GetSolverAndCheckMarker(unsigned short iSolver, unsigned short iInst,
                                          unsigned short iMarker = std::numeric_limits<unsigned short>::max()) const {
    if (iMarker < std::numeric_limits<unsigned short>::max() && iMarker > GetNumberMarkers(iInst)) {
      SU2_MPI::Error("Marker index exceeds size.", CURRENT_FUNCTION);
    }
    auto* solver = solver_container[selected_zone][iInst][MESH_0][iSolver];
    if (solver == nullptr) SU2_MPI::Error("The selected solver does not exist.", CURRENT_FUNCTION);
    return solver;
  }
  unsigned long GetNodeGlobalIndex(unsigned long iPoint, unsigned short iInst) const {
    if (iPoint >= GetNumberNodes(iInst)) {
      SU2_MPI::Error("Node index exceeds mesh size.", CURRENT_FUNCTION);
    }
    return geometry_container[ZONE_0][iInst][MESH_0]->nodes->GetGlobalIndex(iPoint);
  }
  unsigned long GetNumberNodes(unsigned short iInst) const {
    return geometry_container[ZONE_0][iInst][MESH_0]->GetnPoint();
  }

  vector<string> GetDeformableMarkerTags(unsigned short iInst) const {
    const auto nMarker = main_config->GetnMarker_Deform_Mesh();
    vector<string> tags(nMarker);

    for (auto iMarker = 0u; iMarker < nMarker; iMarker++) {
      tags[iMarker] = main_config->GetMarker_Deform_Mesh_TagBound(iMarker);
    }
    return tags;
  }
  vector<string> GetMarkerTags(unsigned short iInst) const {
  const auto nMarker = main_config->GetnMarker_All();
  vector<string> tags(nMarker);

  for (auto iMarker = 0u; iMarker < nMarker; iMarker++) {
    tags[iMarker] = main_config->GetMarker_All_TagBound(iMarker);
  }
  return tags;
}

bool GetNodeDomain(unsigned long iPoint,unsigned short iInst) const {
  if (iPoint >= GetNumberNodes(iInst)) {
    SU2_MPI::Error("Node index exceeds mesh size.", CURRENT_FUNCTION);
  }
  return geometry_container[ZONE_0][iInst][MESH_0]->nodes->GetDomain(iPoint);
}

  inline CPyWrapperMatrixView InitialCoordinates(unsigned short iInst) const {
    if (!main_config->GetDeform_Mesh()) {
      SU2_MPI::Error("Initial coordinates are only available with DEFORM_MESH= YES", CURRENT_FUNCTION);
    }
    auto* coords =
        const_cast<su2activematrix*>(solver_container[selected_zone][iInst][MESH_0][MESH_SOL]->GetNodes()->GetMesh_Coord());
    return CPyWrapperMatrixView(*coords, "InitialCoordinates", true);
  }
 unsigned short GetnInst(){return nInstHB;}
};
