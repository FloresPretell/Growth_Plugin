/*
 * Native CONFORMING (red-green, NO hanging nodes) refinement of the band around a
 * level-set interface, with the new vertices projected onto φ=0 by the
 * LevelSetProjector. This is the projector's payoff for the ghost-fluid flux: the
 * flux disc asserts on hanging-node meshes (the multigrid injection breaks on
 * constrained DoFs), so it needs a CONFORMING mesh with a fine, interface-resolved
 * band — exactly what UG4's flat-grid Refine(grid, sel, projector) produces
 * ("refined in a way to avoid hanging nodes"). Replaces the external MMG remesh
 * with an in-process, projector-driven mesh (fast enough to re-run per step for a
 * moving interface).
 *
 * 2D: refines FACES in the band and closure-refines neighbours to stay conforming.
 * File-based (load → refine npasses → save) so the result reloads as a clean
 * single-level Domain for the flux (no MultiGrid hanging-constraint state to carry).
 */
#ifndef __H__UG_level_set_conforming_refine
#define __H__UG_level_set_conforming_refine

#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include "lib_grid/lg_base.h"
#include "lib_grid/file_io/file_io.h"
#include "lib_grid/tools/selector_grid.h"
#include "lib_grid/refinement/regular_refinement.h"
#include "common/util/smart_pointer.h"
#include "level_set_projector.h"

namespace ug {
namespace LevelSet {

/// Conforming, projector-driven band refinement around an analytic sphere/circle
/// φ = |x−c| − R. Loads inFile, refines `npasses` times (each pass: select the band
/// faces |φ_centre|≤band, red-green Refine with a LevelSetProjector placing new band
/// vertices on φ=0), saves outFile. Returns the vertex count. 2D (faces).
inline int RefineConformingLevelSet(const char* inFile, const char* outFile,
                                    const vector3& center, number R,
                                    number band, int npasses)
{
	Grid grid(GRIDOPT_STANDARD_INTERCONNECTION);
	SubsetHandler sh(grid);
	if (!LoadGridFromFile(grid, sh, inFile))
		UG_THROW("RefineConformingLevelSet: cannot load '" << inFile << "'");

	if (!grid.has_vertex_attachment(aPosition))
		grid.attach_to_vertices(aPosition);
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);

	for (int pass = 0; pass < npasses; ++pass)
	{
	//	projector from analytic φ; project a little beyond the marked band so the
	//	closure vertices near the interface still land on φ=0
		LevelSetProjector proj;
		proj.set_sphere(center, R);
		proj.set_band(band > 0 ? band * 1.5 : -1);
		proj.set_geometry(make_sp(new Geometry<3, 3>(grid, aPosition)));

	//	select the band faces (|φ at the face centre| ≤ band)
		Selector sel(grid);
		for (FaceIterator it = grid.begin<Face>(); it != grid.end<Face>(); ++it)
		{
			Face* f = *it;
			const size_t n = f->num_vertices();
			vector3 c(0, 0, 0);
			for (size_t i = 0; i < n; ++i) VecAdd(c, c, aaPos[f->vertex(i)]);
			VecScale(c, c, 1.0 / (number)n);
			vector3 d; VecSubtract(d, c, center);
			if (std::fabs(VecLength(d) - R) <= band) sel.select(f);
		}

		if (sel.num<Face>() == 0) break;   // nothing left in the band
		Refine(grid, sel, &proj);
	}

	if (!SaveGridToFile(grid, sh, outFile))
		UG_THROW("RefineConformingLevelSet: cannot save '" << outFile << "'");
	return (int)grid.num_vertices();
}

// ---------------------------------------------------------------------------
//  Conforming (no-hanging) refinement of a TUBE/TREE morphology given by an SWC
//  skeleton. phi(x) = min over SWC segments of [dist(x, segment) - r(segment)];
//  each pass GRADED (band halves) so the interior is uniformly fine and the
//  exterior steps down one level -> 2:1 balanced, red-green conforming (no fans).
//  NO projector (a level-set simulation keeps the interface implicit in phi, so
//  the grid stays clean through curves). File-based; result reloads as a clean
//  single-level Domain the ghost-fluid flux can assemble on.
// ---------------------------------------------------------------------------
struct LSSeg { vector3 a, b; number ra, rb; };

inline number LSSegDist(const vector3& x, const LSSeg& s)
{
	vector3 ab; VecSubtract(ab, s.b, s.a);
	vector3 ax; VecSubtract(ax, x, s.a);
	const number dd = VecDot(ab, ab);
	number t = dd > 0 ? VecDot(ax, ab) / dd : 0;
	if (t < 0) t = 0; else if (t > 1) t = 1;
	vector3 c; VecScaleAdd(c, 1.0, s.a, t, ab);
	return VecDistance(x, c) - (s.ra + t * (s.rb - s.ra));
}

inline bool LSLoadSWC(const char* swcFile, std::vector<LSSeg>& seg, bool planar)
{
	std::ifstream f(swcFile);
	if (!f) return false;
	std::vector<vector3> pos; std::vector<number> rad; std::vector<int> par;
	std::string line;
	while (std::getline(f, line)) {
		if (line.empty() || line[0] == '#') continue;
		std::istringstream ss(line);
		int id, ty, p; number x, y, z, r;
		if (!(ss >> id >> ty >> x >> y >> z >> r >> p)) continue;
		if ((int)pos.size() < id + 1) { pos.resize(id + 1); rad.resize(id + 1); par.resize(id + 1, -1); }
		pos[id] = vector3(x, y, planar ? 0.0 : z); rad[id] = r; par[id] = p;
	}
	for (size_t id = 0; id < par.size(); ++id) {
		int p = par[id];
		if (p > 0 && p < (int)pos.size()) {
			LSSeg s; s.a = pos[p]; s.b = pos[id]; s.ra = rad[p]; s.rb = rad[id];
			seg.push_back(s);
		}
	}
	return !seg.empty();
}

/// Conforming SDF refinement of an SWC morphology. Returns the vertex count.
inline int RefineConformingSWC(const char* inFile, const char* swcFile, const char* outFile,
                               number band, int npasses, bool planar)
{
	Grid grid(GRIDOPT_STANDARD_INTERCONNECTION);
	SubsetHandler sh(grid);
	if (!LoadGridFromFile(grid, sh, inFile))
		UG_THROW("RefineConformingSWC: cannot load '" << inFile << "'");
	std::vector<LSSeg> seg;
	if (!LSLoadSWC(swcFile, seg, planar))
		UG_THROW("RefineConformingSWC: cannot load SWC '" << swcFile << "'");

	if (!grid.has_vertex_attachment(aPosition)) grid.attach_to_vertices(aPosition);
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);

	for (int pass = 0; pass < npasses; ++pass)
	{
		const number bk = band * std::pow(0.5, pass);   // graded 2:1
		Selector sel(grid);
		for (FaceIterator it = grid.begin<Face>(); it != grid.end<Face>(); ++it)
		{
			Face* fc = *it;
			const size_t n = fc->num_vertices();
			vector3 c(0, 0, 0);
			for (size_t i = 0; i < n; ++i) VecAdd(c, c, aaPos[fc->vertex(i)]);
			VecScale(c, c, 1.0 / (number)n);
			number phi = 1e30;
			for (size_t k = 0; k < seg.size(); ++k) { number d = LSSegDist(c, seg[k]); if (d < phi) phi = d; }
			if (phi <= bk) sel.select(fc);
		}
		if (sel.num<Face>() == 0) break;
		Refine(grid, sel, NULL);   // red-green conforming, no projector (implicit interface)
	}
	if (!SaveGridToFile(grid, sh, outFile))
		UG_THROW("RefineConformingSWC: cannot save '" << outFile << "'");
	return (int)grid.num_vertices();
}

}  // namespace LevelSet
}  // namespace ug

#endif
