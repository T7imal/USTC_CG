#include <Engine/MeshEdit/MinSurf.h>

#include <Engine/Primitive/TriMesh.h>

#include <Eigen/Sparse>

using namespace Ubpa;

using namespace std;
using namespace Eigen;

MinSurf::MinSurf(Ptr<TriMesh> triMesh)
	: heMesh(make_shared<HEMesh<V>>()) {
	Init(triMesh);
}

void MinSurf::Clear() {
	heMesh->Clear();
	triMesh = nullptr;
}

bool MinSurf::Init(Ptr<TriMesh> triMesh) {
	Clear();

	if (triMesh == nullptr)
		return true;

	if (triMesh->GetType() == TriMesh::INVALID) {
		printf("ERROR::MinSurf::Init:\n"
			"\t""trimesh is invalid\n");
		return false;
	}

	// init half-edge structure
	size_t nV = triMesh->GetPositions().size();
	vector<vector<size_t>> triangles;
	triangles.reserve(triMesh->GetTriangles().size());
	for (auto triangle : triMesh->GetTriangles())
		triangles.push_back({ triangle->idx[0], triangle->idx[1], triangle->idx[2] });
	heMesh->Reserve(nV);
	heMesh->Init(triangles);

	if (!heMesh->IsTriMesh() || !heMesh->HaveBoundary()) {
		printf("ERROR::MinSurf::Init:\n"
			"\t""trimesh is not a triangle mesh or hasn't a boundaries\n");
		heMesh->Clear();
		return false;
	}

	// triangle mesh's positions ->  half-edge structure's positions
	for (int i = 0; i < nV; i++) {
		auto v = heMesh->Vertices().at(i);
		v->pos = triMesh->GetPositions()[i].cast_to<vecf3>();
	}

	this->triMesh = triMesh;
	return true;
}

bool MinSurf::Run() {
	if (heMesh->IsEmpty() || !triMesh) {
		printf("ERROR::MinSurf::Run\n"
			"\t""heMesh->IsEmpty() || !triMesh\n");
		return false;
	}

	Minimize();

	// half-edge structure -> triangle mesh
	size_t nV = heMesh->NumVertices();
	size_t nF = heMesh->NumPolygons();
	vector<pointf3> positions;
	vector<unsigned> indice;
	positions.reserve(nV);
	indice.reserve(3 * nF);
	for (auto v : heMesh->Vertices())
		positions.push_back(v->pos.cast_to<pointf3>());
	for (auto f : heMesh->Polygons()) { // f is triangle
		for (auto v : f->BoundaryVertice()) // vertices of the triangle
			indice.push_back(static_cast<unsigned>(heMesh->Index(v)));
	}

	triMesh->Init(indice, positions);

	return true;
}

void MinSurf::Minimize() {
	// TODO
	size_t nV = heMesh->NumVertices();
	// 计算微分坐标左矩阵
	vector<Triplet<double>> L;
	for (size_t i = 0; i < nV; i++) {
		auto vi = heMesh->Vertices().at(i);
		L.push_back(Triplet<double>(i, i, 1));
		if (!vi->IsBoundary()) {
			for (auto vj : vi->AdjVertices()) {
				L.push_back(Triplet<double>(i, heMesh->Index(vj), -1.0 / vi->Degree()));
			}
		}
	}
	SparseMatrix<double> Lmat(nV, nV);
	Lmat.setFromTriplets(L.begin(), L.end());
	Lmat.makeCompressed();

	SimplicialLDLT<SparseMatrix<double>> solver;
	solver.compute(Lmat);
	if (solver.info() != Success) {
		cout << "WARNING::MinSurf::Minimize:" << endl
			<< "\t" << "solver decompose fail" << endl;
		return;
	}

	// 计算三个维度右矩阵
	VectorXd bx(nV), by(nV), bz(nV);
	bx.setZero(); by.setZero(); bz.setZero();
	VectorXd x(nV), y(nV), z(nV);
	x.setZero(); y.setZero(); z.setZero();

	//固定边界
	for (size_t i = 0; i < nV; i++) {
		auto vi = heMesh->Vertices().at(i);
		if (vi->IsBoundary()) {
			bx(i) = vi->pos[0];
			by(i) = vi->pos[1];
			bz(i) = vi->pos[2];
		}
	}

	// 求解坐标
	x = solver.solve(bx);
	y = solver.solve(by);
	z = solver.solve(bz);

	// 更新坐标
	for (size_t i = 0; i < nV; i++) {
		auto vi = heMesh->Vertices().at(i);
		vi->pos.at(0) = x(i);
		vi->pos.at(1) = y(i);
		vi->pos.at(2) = z(i);
	}

	// cout << "WARNING::MinSurf::Minimize:" << endl
	// 	<< "\t" << "not implemented" << endl;
}
