/*
 * Refinement_LS — level-set-driven refinement / interface-fitted projectors for the
 * moving-boundary neuron growth model. RELOCATED 2026-07-16 out of Dmitry's LevelSet
 * plugin (we use LevelSet but do not modify it) into Nicole's Growth_Plugin.
 *
 * The eight worker headers (level_set_projector.h, _projector_tools.h, _conforming_refine.h,
 * _swc_tools.h, _marker_tools.h, _tube_projector.h, _tube_builder.h, _tube_tools.h) are
 * self-contained (only ugcore/lib_disc + each other) and live in this folder. Their symbols
 * stay in namespace ug::LevelSet; this file wires them into the NeuroGrowth registry.
 *
 * Registration is split exactly as it was in the LevelSet plugin:
 *   RegisterDomainAlgebra<TDomain,TAlgebra>  -> the GridFunction/TDomain functions
 *   RegisterCommon                            -> the (dimension-independent) projector classes
 *                                                + the file-based conforming refiners
 */
#ifndef __H__UG_Growth_Refinement_LS
#define __H__UG_Growth_Refinement_LS

#include <string>
#include "bridge/util.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_grid/refinement/projectors/refinement_projector.h"

#include "level_set_projector.h"
#include "level_set_projector_tools.h"
#include "level_set_marker_tools.h"
#include "level_set_conforming_refine.h"
#include "level_set_tube_projector.h"
#include "level_set_tube_builder.h"
#include "level_set_tube_tools.h"
#include "level_set_swc_tools.h"

namespace ug {
namespace RefinementLS {

//	GridFunction / TDomain dependent functions (was LevelSet::Functionality::DomainAlgebra)
template <typename TDomain, typename TAlgebra>
void RegisterDomainAlgebra(bridge::Registry& reg, std::string grp)
{
	typedef GridFunction<TDomain, TAlgebra> function_type;

//	LevelSetProjector: read a level-set GridFunction φ into the projector (discrete v2)
	reg.add_function("SetLevelSetProjectorLSF",
	                 &LevelSet::SetLevelSetProjectorLSF<function_type>, grp,
	                 "", "projector#levelSetGF#component",
	                 "Fill a LevelSetProjector from a level-set GridFunction");

//	LevelSetTubeProjector: derive per-vertex tube coords for all grid vertices from the SWC
	reg.add_function("PrepareLevelSetTube",
	                 &LevelSet::PrepareLevelSetTube<TDomain>, grp,
	                 "", "projector#domain",
	                 "Derive per-vertex (axial,angular,radial) by projecting onto the SWC centreline");

//	φ-driven markers (orientation of refinement along/across ∇φ; growth band; coarsen)
	reg.add_function("MarkLevelSetTube",
	                 &LevelSet::MarkLevelSetTube<function_type>, grp,
	                 "num marked", "refiner#u#cmp#outerBand",
	                 "Mark the whole interior (phi<=outerBand) for refinement — a solid tube");
	reg.add_function("MarkLevelSetGrowthBand",
	                 &LevelSet::MarkLevelSetGrowthBand<function_type>, grp,
	                 "num marked", "refiner#uOld#uNew#cmp#band#collar",
	                 "Mark only the interior the interface just swept into (growth): "
	                 "phiOld>band AND phiNew<=band+collar — the moving-interface refiner");
	reg.add_function("MarkLevelSetCoarsen",
	                 &LevelSet::MarkLevelSetCoarsen<function_type>, grp,
	                 "num marked", "refiner#u#cmp#wc",
	                 "Mark elements with |phi_centre|>=wc for RM_COARSEN (fast C++ coarsen marker)");
	reg.add_function("InterpolateSWCPhi",
	                 &LevelSet::InterpolateSWCPhi<function_type>, grp,
	                 "", "u#cmp#swcFile#radiusOffset#planar",
	                 "Fill phi = min dist(SWC segments) - radius - radiusOffset in C++ (fast, no Lua loop)");
	reg.add_function("FillSignEps",
	                 &LevelSet::FillSignEps<function_type>, grp,
	                 "", "dst#src#cmp#eps",
	                 "Fill dst = src/sqrt(src^2+eps^2) (smoothed sign) — the frozen reinit source sign(phi0)");
	reg.add_function("AddSignSource",
	                 &LevelSet::AddSignSource<function_type>, grp,
	                 "", "u#phi0#cmp#dtau#eps",
	                 "Operator-split reinit source: u += dtau*sign_eps(phi0) (interface preserved)");
	reg.add_function("MarkAnisotropicLevelSet",
	                 &LevelSet::MarkAnisotropicLevelSet<function_type>, grp,
	                 "numMarked", "refiner#levelSetGF#component#band#acrossInterface#interior#maxLvl",
	                 "Anisotropically mark elements oriented to the level-set normal; "
	                 "interior=true refines the whole interior (radial shells), false = |phi|<=band shell; "
	                 "maxLvl caps refinement depth (-1 = no cap)");
	reg.add_function("CountInterfaceBelowLevel",
	                 &LevelSet::CountInterfaceBelowLevel<function_type>, grp,
	                 "num cut cells below maxLvl", "levelSetGF#component#maxLvl",
	                 "VERIFY adaptivity: count phi=0-cut surface cells whose level < maxLvl "
	                 "(0 = the moving front is fully refined everywhere)");
	reg.add_function("MaxGradDevInBand",
	                 &LevelSet::MaxGradDevInBand<function_type>, grp,
	                 "max ||grad phi|-1|", "levelSetGF#component#band",
	                 "VERIFY reinit: max over band elements of ||grad phi|-1| "
	                 "(0 = phi is a perfect signed distance; large = needs reinit)");
	reg.add_function("MeanGradDevInBand",
	                 &LevelSet::MeanGradDevInBand<function_type>, grp,
	                 "mean ||grad phi|-1|", "levelSetGF#component#band",
	                 "VERIFY reinit (robust): mean over band elements of ||grad phi|-1|");
	reg.add_function("ReinitStepNodal",
	                 &LevelSet::ReinitStepNodal<function_type>, grp,
	                 "", "u#phi0#cmp#dtau#eps#band",
	                 "One explicit nodal eikonal reinit step u += dtau*sign(phi0)*(1-|grad phi|) "
	                 "(correct HJ update; interface preserved). Use dtau<h_fine; call n times");
}

//	dimension-independent classes + file-based conforming refiners (registered once)
inline void RegisterCommon(bridge::Registry& reg, std::string grp)
{
//	LevelSetProjector — dimension-independent (position-based). Registered once.
	{
		typedef LevelSet::LevelSetProjector T;
		reg.add_class_<T, RefinementProjector>("LevelSetProjector", grp)
			.add_constructor()
			.add_method("set_sphere", &T::set_sphere, "", "center#radius",
			            "analytic phi = |x-c|-R (reproduces SphereProjector)")
			.add_method("set_cylinder", &T::set_cylinder, "", "axisPoint#axisDir#radius",
			            "analytic phi = dist(x,axis)-R")
			.add_method("set_band", &T::set_band, "", "band",
			            "only project |phi_parent|<=band; far field stays linear")
			.add_method("set_newton_iterations", &T::set_newton_iterations, "", "n")
			.add_method("set_shell_spacing", &T::set_shell_spacing, "", "dPhi",
			            "SHELL mode: project interior onto nested phi=k*dPhi shells (0=off)")
			.add_method("set_num_shells", &T::set_num_shells, "", "nShells#R",
			            "SHELL mode: nShells fractional-radius shells across radius R (dPhi=R/n)")
			.add_method("set_max_radial_disp", &T::set_max_radial_disp, "", "frac",
			            "inversion guard: cap vertex move to frac*(local spacing)")
			.set_construct_as_smart_pointer(true);
	}

//	native CONFORMING (no-hanging) projector-driven band refinement — file based
	reg.add_function("RefineConformingLevelSet",
	                 &LevelSet::RefineConformingLevelSet, grp,
	                 "num vertices", "inFile#outFile#center#radius#band#npasses",
	                 "Red-green conforming band refinement around phi=|x-c|-R, new "
	                 "band vertices projected onto phi=0. Flux-runnable (no hanging).");
//	conforming SDF refinement of an SWC tube/tree morphology (flux-runnable)
	reg.add_function("RefineConformingSWC",
	                 &LevelSet::RefineConformingSWC, grp,
	                 "num vertices", "inFile#swcFile#outFile#band#npasses#planar",
	                 "Graded red-green conforming refinement of an SWC morphology "
	                 "(phi=min dist to branches). No hanging nodes -> flux-runnable.");

//	LevelSetTubeProjector — keeps a tube mesh tube-structured under refinement
	{
		typedef LevelSet::LevelSetTubeProjector T;
		reg.add_class_<T, RefinementProjector>("LevelSetTubeProjector", grp)
			.add_constructor()
			.add_method("set_swc", &T::set_swc, "", "swcFile",
			            "build the centreline model from the SWC longest path")
			.add_method("set_swc_ndense", &T::set_swc_ndense, "", "swcFile#ndense",
			            "as set_swc, with an explicit dense-station count")
			.add_method("set_radial_gate", &T::set_radial_gate, "", "rgate",
			            "deform only vertices with fractional radius <= rgate (interior); box exterior kept")
			.add_method("set_max_radial_disp", &T::set_max_radial_disp, "", "frac",
			            "inversion guard: cap deformation to frac*(parent scale)")
			.add_method("set_deform_base", &T::set_deform_base, "", "bool",
			            "deform the base interior onto the tube in prepare() (default true)")
			.set_construct_as_smart_pointer(true);
	}
//	build the tube-topology base mesh from an SWC (C++ replacement of make_shell_tube.py)
	reg.add_function("CreateLevelSetTube",
	                 &LevelSet::CreateLevelSetTube, grp,
	                 "num vertices", "outFile#swcFile#nrad#nang#naxial",
	                 "Sweep a concentric-radial-shell disk along the SWC longest path "
	                 "into prisms; positions from LevelSetTubeProjector. Saves .ugx.");
}

}  // namespace RefinementLS
}  // namespace ug

#endif
