#include "PolynomialList.h"
#include <iostream>
#include <fstream>

using namespace std;

PolynomialList::PolynomialList(const PolynomialList& other) {
	// TODO
	m_Polynomial = other.m_Polynomial;
}

PolynomialList::PolynomialList(const string& file) {
	// TODO
	ReadFromFile(file);
}

PolynomialList::PolynomialList(const double* cof, const int* deg, int n) {
	// TODO
	for (int i = 0; i < n; ++i) {
		m_Polynomial.push_back(Term(deg[i], cof[i]));
	}
}

PolynomialList::PolynomialList(const vector<int>& deg, const vector<double>& cof) {
	// TODO
	int n = min(deg.size(), cof.size());
	for (int i = 0; i < n; ++i) {
		m_Polynomial.push_back(Term(deg[i], cof[i]));
	}
}

double PolynomialList::coff(int i) const {
	// TODO
	for (auto t : m_Polynomial) {
		if (t.deg == i) {
			return t.cof;
		}
	}
	return 0;
}

double& PolynomialList::coff(int i) {
	// TODO
	return AddOneTerm(Term(i, 0)).cof;
}

void PolynomialList::compress() {
	// TODO
	auto itr = m_Polynomial.begin();
	while (itr != m_Polynomial.end()) {
		if (fabs((*itr).cof) < 1.0e-10)
			itr = m_Polynomial.erase(itr);
		else
			itr++;
	}
}

PolynomialList PolynomialList::operator+(const PolynomialList& right) const {
	// TODO
	PolynomialList a = PolynomialList(right);
	for (auto t : this->m_Polynomial) {
		a.m_Polynomial.push_back(t);
	}
	a.compress();
	return a;
}

PolynomialList PolynomialList::operator-(const PolynomialList& right) const {
	// TODO
	PolynomialList a = PolynomialList(right);
	for (auto t : this->m_Polynomial) {
		t.cof = -t.cof;
		a.m_Polynomial.push_back(t);
	}
	a.compress();
	return a;
}

PolynomialList PolynomialList::operator*(const PolynomialList& right) const {
	// TODO
	PolynomialList a = PolynomialList();
	for (auto t : this->m_Polynomial) {
		for (auto x : right.m_Polynomial) {
			double cof = t.cof * x.cof;
			int deg = t.deg + x.deg;
			a.AddOneTerm(Term(deg, cof));
		}
	}
	a.compress();
	return a;
}

PolynomialList& PolynomialList::operator=(const PolynomialList& right) {
	// TODO
	m_Polynomial = right.m_Polynomial;
	return *this;
}

void PolynomialList::Print() const {
	auto itr = m_Polynomial.begin();
	while (itr != m_Polynomial.end()) {
		if (itr != m_Polynomial.begin()) {
			if (itr->cof > 0)
				cout << "+";
		}
		if (!abs(itr->cof - 1) < 1e-10) {
			cout << itr->cof;
		}
		if (itr->deg > 0)
			cout << "x^" << itr->deg;
		itr++;
	}
	cout << endl;
}

bool PolynomialList::ReadFromFile(const string& file) {
	// TODO
	ifstream fin;
	fin.open(file.c_str(), ios::in);
	if (!fin.is_open())return false;
	char c;
	int n;
	fin >> c >> n;
	if (c != 'P')return false;
	for (int i = 0; i < n; ++i) {
		int deg;
		double cof;
		fin >> deg >> cof;
		struct Term x(deg, cof);
		AddOneTerm(x);
	}
	fin.close();
	return true;
}

PolynomialList::Term& PolynomialList::AddOneTerm(const Term& term) {
	// TODO
	m_Polynomial.sort();
	auto itr = m_Polynomial.begin();
	for (; itr != m_Polynomial.end(); itr++) {
		if (itr->deg == term.deg) {
			itr->cof += term.cof;
			return *itr;
		}

		if (itr->deg > term.deg)
			break;
	}
	return *m_Polynomial.insert(itr, term);
}
