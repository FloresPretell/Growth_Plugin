/*
 * LevelSetTubeProjector — a RefinementProjector that keeps a tube/neurite mesh
 * TUBE-STRUCTURED under refinement, the level-set/SWC analogue of the UG4
 * NeuriteProjector. Every vertex carries tube coordinates (axial t, angular φ,
 * radial r) derived from the SWC centreline (r fractional: 1 = surface, 0 =
 * centreline → concentric interior shells). On refinement each new vertex averages
 * its parents' tube coords, is reconstructed on the true tube
 *     pos = C(t) + r·R(t)·(cos φ · N(t) + sin φ · B(t))
 * and the averaged coords are written back onto it — so every level stays a clean
 * concentric-shell tube grid (no box-octree spokes).
 *
 * The "saved function" is the SWC: set_swc() builds a dense centreline model
 * (points C, parallel-transport frame N,B, taper radius R, arc-length param t).
 * prepare() derives each existing vertex's tube coords by projecting its position
 * onto that centreline, so no per-vertex attachment needs to be serialized — the
 * SWC is re-applied on every load.
 */
#ifndef __H__UG_level_set_tube_projector
#define __H__UG_level_set_tube_projector

#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "lib_grid/refinement/projectors/refinement_projector.h"
#include "lib_grid/grid/grid.h"
#include "lib_grid/multi_grid.h"
#include "common/math/ugmath.h"

namespace ug {
namespace LevelSet {

struct TubeParams { number t; number ang; number rad; int branch; };

class LevelSetTubeProjector : public RefinementProjector
{
	public:
		LevelSetTubeProjector() : m_have(false) {}
		LevelSetTubeProjector(SPIGeometry3d geometry)
			: RefinementProjector(geometry), m_have(false) {}

		virtual ~LevelSetTubeProjector() {}

	// ---- build the centreline model from the SWC (longest root->tip path) -------
		void set_swc(const char* swcFile) { set_swc_ndense(swcFile, 400); }
		void set_swc_ndense(const char* swcFile, int ndense)
		{
			std::vector<vector3> pts; std::vector<number> rad;
			if (!load_longest_path(swcFile, pts, rad))
				UG_THROW("LevelSetTubeProjector: cannot read a path from '" << swcFile << "'");
			resample(pts, rad, ndense);       // fills m_C, m_R, m_t (arc-length)
			build_frame();                    // fills m_N, m_B (parallel transport)
			m_have = true;
		}

	//	INTERIOR GATE: only vertices with fractional radius rad <= rgate (i.e. inside
	//	the tube + the surface, φ<=0) are deformed onto the tube; the exterior (box)
	//	is left untouched so the box is kept for the moving level-set interface.
		void set_radial_gate(number rgate)   { m_rgate = rgate; }
		void set_max_radial_disp(number frac) { m_maxRadFrac = frac; }
	//	deform the base interior onto the tube (prepare must run first)
		void set_deform_base(bool b)          { m_deformBase = b; }

	// ---- attach the per-vertex tube-param storage + derive params for all verts --
	//	BOX+SDF mode: derive tube coords for every vertex; DEFORM the interior
	//	(rad<=rgate) onto the tube shells, leave the exterior on the box grid.
		void prepare(Grid& g)
		{
			UG_COND_THROW(!m_have, "LevelSetTubeProjector: set_swc() first.");
			if (!g.has_vertex_attachment(m_aTube)) g.attach_to_vertices(m_aTube);
			m_aaTube.access(g, m_aTube);
			if (!g.has_vertex_attachment(aPosition)) g.attach_to_vertices(aPosition);
			Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);
			for (VertexIterator it = g.begin<Vertex>(); it != g.end<Vertex>(); ++it) {
				const TubeParams p = project_pos(aaPos[*it]);
				m_aaTube[*it] = p;
				if (m_deformBase && p.rad <= m_rgate)     // interior: onto the tube
					aaPos[*it] = reconstruct(p);          // exterior: stay on the box
			}
		}

	// ---- forward map: tube coords -> 3D position --------------------------------
		vector3 reconstruct(const TubeParams& p) const
		{
			number t = p.t; if (t < 0) t = 0; else if (t > 1) t = 1;
			size_t i; number u; locate(t, i, u);
			vector3 C, N, B; number R;
			interp(i, u, C, N, B, R);
			vector3 pos; VecScaleAdd(pos, 1.0, C,
				p.rad * R * std::cos(p.ang), N, p.rad * R * std::sin(p.ang), B);
			return pos;
		}

	// ---- RefinementProjector contract -------------------------------------------
		virtual number new_vertex(Vertex* vrt, Vertex* parent)
		{	m_aaTube[vrt] = m_aaTube[parent]; set_pos(vrt, pos(parent)); return 1; }
		virtual number new_vertex(Vertex* vrt, Edge* parent)   { return place(vrt, parent); }
		virtual number new_vertex(Vertex* vrt, Face* parent)   { return place(vrt, parent); }
		virtual number new_vertex(Vertex* vrt, Volume* parent) { return place(vrt, parent); }

	// ---- one-shot: place an existing vertex from its stored params --------------
		void project(Vertex* vrt) { set_pos(vrt, reconstruct(m_aaTube[vrt])); }

	private:
		template <class TElem>
		vector3 linear_center(TElem* parent, size_t n) const
		{
			vector3 c(0, 0, 0);
			for (size_t i = 0; i < n; ++i) VecAdd(c, c, pos(parent->vertex(i)));
			VecScale(c, c, 1.0 / (number)n);
			return c;
		}

		template <class TElem>
		number place(Vertex* vrt, TElem* parent)
		{
			const size_t n = parent->num_vertices();
			TubeParams avg = average_params(parent, n);
			m_aaTube[vrt] = avg;
			const vector3 xlin = linear_center(parent, n);
			if (avg.rad > m_rgate) {                 // exterior: keep the box midpoint
				set_pos(vrt, xlin);
				return 1;
			}
			vector3 xt = reconstruct(avg);           // interior/surface: onto the tube
			if (m_maxRadFrac > 0) {                   // inversion guard vs the parent scale
				number scale = 0;
				for (size_t i = 0; i < n; ++i) {
					const number d = VecDistance(pos(parent->vertex(i)), xlin);
					if (d > scale) scale = d;
				}
				if (VecDistance(xt, xlin) > m_maxRadFrac * scale) xt = xlin;
			}
			set_pos(vrt, xt);
			return 1;
		}

	//	average parent tube coords: axial mean, angular vector-sum (weighted by
	//	radial), radial = 0.5(min+max), centre-snap when the in-plane vector is
	//	short — mirrors NeuriteProjector::average_params.
		template <class TElem>
		TubeParams average_params(TElem* parent, size_t n) const
		{
			number t = 0, vx = 0, vy = 0, rmin = 1e30, rmax = -1e30;
			int branch = 0;
			for (size_t i = 0; i < n; ++i) {
				const TubeParams& p = m_aaTube[parent->vertex(i)];
				t += p.t;
				vx += p.rad * std::cos(p.ang);
				vy += p.rad * std::sin(p.ang);
				if (p.rad < rmin) rmin = p.rad;
				if (p.rad > rmax) rmax = p.rad;
				branch = p.branch;
			}
			TubeParams out;
			out.t = t / (number)n;
			out.rad = 0.5 * (rmin + rmax);
			out.branch = branch;
			const number vlen = std::sqrt(vx*vx + vy*vy);
			if (vlen / (number)n < 0.4 * (out.rad + 1e-12) && (n == 4 || n == 8)) {
				out.ang = 0.0; out.rad = 0.0;         // snap to the centreline
			} else if (vlen < 1e-8) {
				out.ang = 0.0; out.rad = 0.0;
			} else {
				out.ang = std::atan2(vy, vx);
			}
			return out;
		}

	// ---- inverse map: 3D position -> nearest tube coords (foot-point search) -----
		TubeParams project_pos(const vector3& p) const
		{
			number bestD = 1e30; size_t bi = 0; number bu = 0;
			for (size_t i = 0; i + 1 < m_C.size(); ++i) {
				vector3 seg; VecSubtract(seg, m_C[i+1], m_C[i]);
				const number dd = VecDot(seg, seg);
				vector3 ap; VecSubtract(ap, p, m_C[i]);
				number u = dd > 0 ? VecDot(ap, seg) / dd : 0;
				if (u < 0) u = 0; else if (u > 1) u = 1;
				vector3 foot; VecScaleAdd(foot, 1.0, m_C[i], u, seg);
				const number d = VecDistanceSq(p, foot);
				if (d < bestD) { bestD = d; bi = i; bu = u; }
			}
			vector3 C, N, B; number R;
			interp(bi, bu, C, N, B, R);
			vector3 off; VecSubtract(off, p, C);
			vector3 T; VecSubtract(T, m_C[bi+1], m_C[bi]);
			const number Tl = VecLength(T);
			if (Tl > 1e-30) { VecScale(T, T, 1.0 / Tl); number a = VecDot(off, T);
			                  vector3 ax; VecScale(ax, T, a); VecSubtract(off, off, ax); }
			TubeParams q;
			q.branch = 0;
			q.t = m_t[bi] + bu * (m_t[bi+1] - m_t[bi]);
			const number rlen = VecLength(off);
			q.rad = (R > 1e-30) ? rlen / R : 0.0;
			q.ang = (rlen > 1e-12) ? std::atan2(VecDot(off, B), VecDot(off, N)) : 0.0;
			return q;
		}

	// ---- centreline model helpers ------------------------------------------------
		void locate(number t, size_t& i, number& u) const
		{
			// find segment [i,i+1] with m_t[i] <= t <= m_t[i+1]
			size_t lo = 0, hi = m_t.size() - 1;
			if (t <= m_t.front()) { i = 0; u = 0; return; }
			if (t >= m_t.back())  { i = m_t.size() - 2; u = 1; return; }
			while (hi - lo > 1) { size_t mid = (lo + hi) / 2; (m_t[mid] <= t ? lo : hi) = mid; }
			i = lo;
			const number denom = m_t[i+1] - m_t[i];
			u = denom > 1e-30 ? (t - m_t[i]) / denom : 0.0;
		}
		void interp(size_t i, number u, vector3& C, vector3& N, vector3& B, number& R) const
		{
			VecScaleAdd(C, 1.0 - u, m_C[i], u, m_C[i+1]);
			VecScaleAdd(N, 1.0 - u, m_N[i], u, m_N[i+1]);
			VecScaleAdd(B, 1.0 - u, m_B[i], u, m_B[i+1]);
			R = (1.0 - u) * m_R[i] + u * m_R[i+1];
			number nl = VecLength(N); if (nl > 1e-30) VecScale(N, N, 1.0 / nl);
			// re-orthogonalise B against N
			number d = VecDot(B, N); vector3 corr; VecScale(corr, N, d); VecSubtract(B, B, corr);
			number bl = VecLength(B); if (bl > 1e-30) VecScale(B, B, 1.0 / bl);
		}

	// ---- SWC reading + resample + parallel-transport frame ----------------------
		bool load_longest_path(const char* swcFile, std::vector<vector3>& pts,
		                       std::vector<number>& rad) const
		{
			std::ifstream f(swcFile); if (!f) return false;
			std::vector<vector3> pos; std::vector<number> rr; std::vector<int> par;
			std::string line;
			while (std::getline(f, line)) {
				if (line.empty() || line[0] == '#') continue;
				std::istringstream ss(line);
				int id, ty, p; number x, y, z, r;
				if (!(ss >> id >> ty >> x >> y >> z >> r >> p)) continue;
				if ((int)pos.size() < id + 1) { pos.resize(id+1); rr.resize(id+1); par.resize(id+1, -2); }
				pos[id] = vector3(x, y, z); rr[id] = r; par[id] = p;
			}
			// children
			std::vector<std::vector<int> > kids(pos.size());
			std::vector<int> roots;
			for (size_t id = 0; id < par.size(); ++id) {
				if (par[id] == -2) continue;                     // absent id
				if (par[id] < 0 || (size_t)par[id] >= pos.size()) roots.push_back((int)id);
				else kids[par[id]].push_back((int)id);
			}
			// longest path by arc length (iterative DFS)
			number bestLen = -1; std::vector<int> best;
			std::vector<int> stackId; std::vector<number> stackLen;
			std::vector<std::vector<int> > stackPath;
			for (size_t r = 0; r < roots.size(); ++r) {
				stackId.push_back(roots[r]); stackLen.push_back(0.0);
				stackPath.push_back(std::vector<int>(1, roots[r]));
			}
			while (!stackId.empty()) {
				int id = stackId.back(); stackId.pop_back();
				number len = stackLen.back(); stackLen.pop_back();
				std::vector<int> path = stackPath.back(); stackPath.pop_back();
				if (kids[id].empty()) {
					if (len > bestLen) { bestLen = len; best = path; }
					continue;
				}
				for (size_t c = 0; c < kids[id].size(); ++c) {
					int ch = kids[id][c];
					number seg = VecDistance(pos[id], pos[ch]);
					std::vector<int> np = path; np.push_back(ch);
					stackId.push_back(ch); stackLen.push_back(len + seg); stackPath.push_back(np);
				}
			}
			if (best.size() < 2) return false;
			for (size_t k = 0; k < best.size(); ++k) { pts.push_back(pos[best[k]]); rad.push_back(rr[best[k]]); }
			return true;
		}

		void resample(const std::vector<vector3>& pts, const std::vector<number>& rad, int ndense)
		{
			std::vector<number> s(pts.size(), 0.0);
			for (size_t i = 1; i < pts.size(); ++i) s[i] = s[i-1] + VecDistance(pts[i-1], pts[i]);
			const number total = s.back();
			m_C.clear(); m_R.clear(); m_t.clear();
			if (total <= 0) { m_C = pts; m_R = rad; m_t.assign(pts.size(), 0.0); return; }
			for (int k = 0; k < ndense; ++k) {
				const number su = total * (number)k / (number)(ndense - 1);
				size_t j = 0; while (j + 1 < s.size() && s[j+1] < su) ++j;
				const number denom = s[j+1] - s[j];
				const number u = denom > 1e-30 ? (su - s[j]) / denom : 0.0;
				vector3 C; VecScaleAdd(C, 1.0 - u, pts[j], u, pts[j+1]);
				m_C.push_back(C);
				m_R.push_back((1.0 - u) * rad[j] + u * rad[j+1]);
				m_t.push_back(su / total);
			}
		}

		void build_frame()
		{
			const size_t n = m_C.size();
			m_N.assign(n, vector3(0,0,0)); m_B.assign(n, vector3(0,0,0));
			std::vector<vector3> T(n);
			for (size_t i = 0; i < n; ++i) {
				vector3 t;
				if (i == 0)          VecSubtract(t, m_C[1], m_C[0]);
				else if (i == n-1)   VecSubtract(t, m_C[n-1], m_C[n-2]);
				else                 VecSubtract(t, m_C[i+1], m_C[i-1]);
				number l = VecLength(t); if (l > 1e-30) VecScale(t, t, 1.0/l);
				T[i] = t;
			}
			// seed normal orthogonal to T[0]
			vector3 seed = (std::fabs(T[0].z()) < 0.9) ? vector3(0,0,1) : vector3(1,0,0);
			number d = VecDot(seed, T[0]); vector3 nn; VecScaleAdd(nn, 1.0, seed, -d, T[0]);
			number nl = VecLength(nn); if (nl > 1e-30) VecScale(nn, nn, 1.0/nl);
			m_N[0] = nn;
			for (size_t i = 1; i < n; ++i) {
				// rotate previous normal by the frame rotation T[i-1]->T[i]
				vector3 v; VecCross(v, T[i-1], T[i]);
				number s = VecLength(v), c = VecDot(T[i-1], T[i]);
				vector3 ni = m_N[i-1];
				if (s > 1e-9) {
					VecScale(v, v, 1.0/s); number ang = std::atan2(s, c);
					vector3 term1; VecScale(term1, ni, std::cos(ang));
					vector3 vxn; VecCross(vxn, v, ni); VecScale(vxn, vxn, std::sin(ang));
					vector3 term3; VecScale(term3, v, VecDot(v, ni) * (1.0 - std::cos(ang)));
					VecAdd(ni, term1, vxn); VecAdd(ni, ni, term3);
				}
				number dd = VecDot(ni, T[i]); vector3 corr; VecScale(corr, T[i], dd); VecSubtract(ni, ni, corr);
				number l = VecLength(ni); if (l > 1e-30) VecScale(ni, ni, 1.0/l);
				m_N[i] = ni;
			}
			for (size_t i = 0; i < n; ++i) { vector3 b; VecCross(b, T[i], m_N[i]); m_B[i] = b; }
		}

	// ---- data -------------------------------------------------------------------
		std::vector<vector3> m_C, m_N, m_B;   // dense centreline points + PT frame
		std::vector<number>  m_R, m_t;        // taper radius, normalised arc-length param
		bool m_have;
		number m_rgate      = 1.05;           // deform only rad<=rgate (interior+surface)
		number m_maxRadFrac = 0.5;            // inversion guard (fraction of parent scale)
		bool   m_deformBase = true;           // deform the base interior in prepare()
		Attachment<TubeParams> m_aTube;
		Grid::VertexAttachmentAccessor<Attachment<TubeParams> > m_aaTube;
};

}  // namespace LevelSet
}  // namespace ug

#endif
