/*
 * Domain bridge for the LevelSetTubeProjector: derive per-vertex tube coordinates
 * for every vertex of the domain's grid by projecting onto the SWC centreline. Call
 * after loading the tube .ugx and before attaching the projector to the refiner.
 */
#ifndef __H__UG_level_set_tube_tools
#define __H__UG_level_set_tube_tools

#include "level_set_tube_projector.h"

namespace ug {
namespace LevelSet {

template <typename TDomain>
void PrepareLevelSetTube(SmartPtr<LevelSetTubeProjector> proj, SmartPtr<TDomain> dom)
{
	proj->prepare(*dom->grid());
}

}  // namespace LevelSet
}  // namespace ug

#endif
