#pragma once

#include <Basic/HeapObj.h>

namespace Ubpa {
	class TriMesh;
	class MinSurf;

	// mesh boundary == 1
	class Paramaterize : public HeapObj {
	public:
		Paramaterize(Ptr<TriMesh> triMesh);
	public:
		static const Ptr<Paramaterize> New(Ptr<TriMesh> triMesh) {
			return Ubpa::New<Paramaterize>(triMesh);
		}
	public:
		enum Boundary_Type {
			kCircle,
			kSquare
		}boundary_type;

		enum Barycentric_Type {
			kUniform,
			kCotangent
		}barycentric_type;

		enum Display_Status {
			koff,
			kon
		}display_status;

		void Set_Method(Barycentric_Type barycentric);
		void Set_Boundary(Boundary_Type boundary);
		void Parameterize();
	public:
		void Clear();
		bool Init(Ptr<TriMesh> triMesh);

		bool Run();
	private:
		class V;
		class E;
		class P;
		class V : public TVertex<V, E, P> {
		public:
			vecf3 pos;
		};
		class E : public TEdge<V, E, P> {};
		class P :public TPolygon<V, E, P> {};
	private:
		Ptr<TriMesh> triMesh;
		const Ptr<HEMesh<V>> heMesh; // vertice order is same with triMesh
	};
}
