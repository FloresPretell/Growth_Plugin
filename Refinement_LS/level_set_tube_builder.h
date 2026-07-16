/*
 * CreateLevelSetTube — the C++ builder of the tube-TOPOLOGY base mesh from an SWC,
 * the proper-UG4 replacement of make_shell_tube.py. A triangulated disk cross-
 * section (centre + nrad concentric rings × nang sectors) is swept along the SWC
 * centreline into PRISMS, giving concentric interior shells by construction. Vertex
 * positions come from the LevelSetTubeProjector's forward map (one source of truth),
 * so the base mesh and every later refinement level share the same tube geometry.
 * The saved .ugx reloads as a clean Domain; the driver re-attaches a
 * LevelSetTubeProjector (built from the same SWC) and refines.
 */
#ifndef __H__UG_level_set_tube_builder
#define __H__UG_level_set_tube_builder

#include <vector>
#include <cmath>
#include "lib_grid/lg_base.h"
#include "lib_grid/file_io/file_io.h"
#include "lib_grid/grid_objects/grid_objects.h"
#include "level_set_tube_projector.h"

namespace ug {
namespace LevelSet {

/// Build a tube-topology shell mesh from an SWC (longest path) and save it.
/// nrad = concentric radial layers, nang = angular sectors, naxial = axial steps.
inline int CreateLevelSetTube(const char* outFile, const char* swcFile,
                              int nrad, int nang, int naxial)
{
	LevelSetTubeProjector proj;
	proj.set_swc(swcFile);

	Grid grid(GRIDOPT_STANDARD_INTERCONNECTION | GRIDOPT_AUTOGENERATE_SIDES);
	SubsetHandler sh(grid);
	grid.attach_to_vertices(aPosition);
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);

	const int per = 1 + nrad * nang;                 // nodes per cross-section
	std::vector<Vertex*> V;
	V.reserve((naxial + 1) * per);

	// ---- vertices: centre + concentric rings, placed via the projector ----------
	for (int s = 0; s <= naxial; ++s) {
		const number t = (number)s / (number)naxial;
		TubeParams c; c.t = t; c.ang = 0; c.rad = 0; c.branch = 0;
		Vertex* vc = *grid.create<RegularVertex>(); aaPos[vc] = proj.reconstruct(c); V.push_back(vc);
		for (int k = 1; k <= nrad; ++k) {
			for (int j = 0; j < nang; ++j) {
				TubeParams p; p.t = t; p.ang = 2.0 * PI * (number)j / (number)nang;
				p.rad = (number)k / (number)nrad; p.branch = 0;
				Vertex* v = *grid.create<RegularVertex>(); aaPos[v] = proj.reconstruct(p); V.push_back(v);
			}
		}
	}
	// section-local node index -> global vertex
	#define NID(s, local) V[(s) * per + (local)]
	// disk-local index of ring k (1..nrad), sector j
	// (k==0 -> centre, local 0)
	struct L { static int idx(int k, int j, int nang) { return k == 0 ? 0 : 1 + (k-1)*nang + (j % nang); } };

	// ---- disk triangulation (centre fan + radial strips) ------------------------
	std::vector<int> triA, triB, triC;   // disk-local indices
	for (int j = 0; j < nang; ++j) {
		triA.push_back(L::idx(0, 0, nang));
		triB.push_back(L::idx(1, j, nang));
		triC.push_back(L::idx(1, j+1, nang));
	}
	for (int k = 1; k < nrad; ++k) {
		for (int j = 0; j < nang; ++j) {
			int a0 = L::idx(k, j, nang),   a1 = L::idx(k, j+1, nang);
			int b0 = L::idx(k+1, j, nang), b1 = L::idx(k+1, j+1, nang);
			triA.push_back(a0); triB.push_back(b0); triC.push_back(b1);
			triA.push_back(a0); triB.push_back(b1); triC.push_back(a1);
		}
	}

	// ---- sweep each disk triangle into a prism between s and s+1 -----------------
	for (int s = 0; s < naxial; ++s) {
		for (size_t tIdx = 0; tIdx < triA.size(); ++tIdx) {
			Vertex* v0 = NID(s,   triA[tIdx]); Vertex* v1 = NID(s,   triB[tIdx]); Vertex* v2 = NID(s,   triC[tIdx]);
			Vertex* v3 = NID(s+1, triA[tIdx]); Vertex* v4 = NID(s+1, triB[tIdx]); Vertex* v5 = NID(s+1, triC[tIdx]);
			grid.create<Prism>(PrismDescriptor(v0, v1, v2, v3, v4, v5));
		}
	}
	#undef NID

	// ---- subsets: everything -> "Inner"; "Boundary" declared (empty) ------------
	for (VertexIterator it = grid.begin<Vertex>(); it != grid.end<Vertex>(); ++it) sh.assign_subset(*it, 0);
	for (EdgeIterator   it = grid.begin<Edge>();   it != grid.end<Edge>();   ++it) sh.assign_subset(*it, 0);
	for (FaceIterator   it = grid.begin<Face>();   it != grid.end<Face>();   ++it) sh.assign_subset(*it, 0);
	for (VolumeIterator it = grid.begin<Volume>(); it != grid.end<Volume>(); ++it) sh.assign_subset(*it, 0);
	sh.subset_info(0).name = "Inner";
	sh.subset_info(1).name = "Boundary";

	if (!SaveGridToFile(grid, sh, outFile))
		UG_THROW("CreateLevelSetTube: cannot save '" << outFile << "'");
	return (int)grid.num_vertices();
}

}  // namespace LevelSet
}  // namespace ug

#endif
