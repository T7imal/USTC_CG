#include <Engine/MeshEdit/Paramaterize.h>

#include <Engine/MeshEdit/MinSurf.h>

#include <Engine/Primitive/TriMesh.h>
#include <Eigen/Sparse>

using namespace Ubpa;
using namespace Eigen;
using namespace std;

#define PI 3.14159265358979323846

Paramaterize::Paramaterize(Ptr<TriMesh> triMesh) {
	// TODO
	Init(triMesh);

}

void Paramaterize::Clear() {
	// TODO
	heMesh->Clear();
	triMesh = nullptr;

}

bool Paramaterize::Init(Ptr<TriMesh> triMesh) {
	// TODO
	Clear();

	if (triMesh == nullptr)
		return true;

	if (triMesh->GetType() == TriMesh::INVALID) {
		printf("ERROR::Paramaterize::Init:\n"
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
		printf("ERROR::Paramaterize::Init:\n"
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


	return false;
}

bool Paramaterize::Run() {
	// TODO
	if (heMesh->IsEmpty() || !triMesh) {
		printf("ERROR::Paramaterize::Run\n"
			"\t""heMesh->IsEmpty() || !triMesh\n");
		return false;
	}

	Set_Boundary(kCircle);
	Set_Method(kUniform);
	Parameterize();

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

void Paramaterize::Set_Method(Barycentric_Type barycentric) {
	barycentric_type = barycentric;
}

void Paramaterize::Set_Boundary(Boundary_Type boundary) {
	boundary_type = boundary;
}

void Paramaterize::Parameterize() {
	// 计算左矩阵
	vector<Triplet<double>> T;
	size_t nV = heMesh->NumVertices();

	switch (barycentric_type) {
	case kUniform:
		for (size_t i = 0; i < nV; i++) {
			auto vi = heMesh->Vertices().at(i);
			T.push_back(Triplet<double>(i, i, 1));
			if (!vi->IsBoundary()) {
				for (auto vj : vi->AdjVertices()) {
					T.push_back(Triplet<double>(i, heMesh->Index(vj), -1.0 / vi->Degree()));
				}
			}
		}
		break;
	case kCotangent:
		for (size_t i = 0; i < nV; i++) {
			auto vi = heMesh->Vertices().at(i);
			T.push_back(Triplet<double>(i, i, 1));
			if (!vi->IsBoundary()) {
				double sum = 0;
				vector<double> weights;
				for (size_t j = 0; j < vi->Degree(); j++) {
					auto vj = vi->AdjVertices()[j];
					auto vprev = vi->AdjVertices()[(j + vi->Degree() - 1) % vi->Degree()];
					auto vnext = vi->AdjVertices()[(j + 1) % vi->Degree()];
					double cos_alpha = (vj->pos - vi->pos).cos_theta(vprev->pos - vi->pos);
					double cos_beta = (vj->pos - vi->pos).cos_theta(vnext->pos - vi->pos);
					double cot_alpha = cos_alpha / sqrt(1 - cos_alpha * cos_alpha);
					double cot_beta = cos_beta / sqrt(1 - cos_beta * cos_beta);
					sum += cot_alpha + cot_beta;
					weights.push_back(cot_alpha + cot_beta);
				}
				for (size_t j = 0; j < vi->Degree(); j++) {
					auto vj = vi->AdjVertices()[j];
					T.push_back(Triplet<double>(i, heMesh->Index(vj), -weights[j] / sum));
				}
			}
		}
	}

	SparseMatrix<double> L(nV, nV);
	L.setFromTriplets(T.begin(), T.end());
	L.makeCompressed();

	SimplicialLDLT<SparseMatrix<double>> solver;
	solver.compute(L);

	// 计算右矩阵
	VectorXd bx(nV), by(nV), x(nV), y(nV);
	bx.setZero(); by.setZero();
	x.setZero(); y.setZero();

	// 固定边界
	int nB = heMesh->Boundaries().size();
	int n = 0;
	switch (boundary_type) {
	case kCircle:
		n = 0;
		for (size_t i = 0; i < nV; i++) {
			auto vi = heMesh->Vertices().at(i);
			if (vi->IsBoundary()) {
				double theta = 2 * PI * n / nB;
				bx(i) = 0.5 * cos(theta) + 0.5;
				by(i) = 0.5 * sin(theta) + 0.5;
				n++;
			}
		}
		break;
	case kSquare:
		n = 0;
		for (size_t i = 0; i < nV; i++) {
			auto vi = heMesh->Vertices().at(i);
			if (vi->IsBoundary()) {
				if (n < nB / 4) {
					bx(i) = 0;
					by(i) = 4.0 * n / nB;
				}
				else if (n < nB / 2) {
					bx(i) = 4.0 * (n - nB / 4) / nB;
					by(i) = 1;
				}
				else if (n < 3 * nB / 4) {
					bx(i) = 1;
					by(i) = 1 - 4.0 * (n - nB / 2) / nB;
				}
				else {
					bx(i) = 1 - 4.0 * (n - 3 * nB / 4) / nB;
					by(i) = 0;
				}
				n++;
			}
		}
		break;
	}
	x = solver.solve(bx);
	y = solver.solve(by);
	for (size_t i = 0; i < nV; i++) {
		auto vi = heMesh->Vertices().at(i);
		vi->pos[0] = x(i);
		vi->pos[1] = y(i);
		vi->pos[2] = 0;
	}
}