#include "PolynomialMap.h"

#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>

using namespace std;

PolynomialMap::PolynomialMap(const PolynomialMap& other) {
	// TODO
	m_Polynomial = other.m_Polynomial;
}

PolynomialMap::PolynomialMap(const string& file) {
	ReadFromFile(file);
}

PolynomialMap::PolynomialMap(const double* cof, const int* deg, int n) {
	// TODO
	for (int i = 0; i < n; ++i) {
		coff(deg[i]) = cof[i];
	}
}

PolynomialMap::PolynomialMap(const vector<int>& deg, const vector<double>& cof) {
	assert(deg.size() == cof.size());
	// TODO
	for (int i = 0; i < deg.size(); ++i) {
		coff(deg[i]) = cof[i];
	}
}

double PolynomialMap::coff(int i) const {
	// TODO
	auto target = m_Polynomial.find(i);
	if (target == m_Polynomial.end()) {
		return 0.;
	}
	return target->second;
}

double& PolynomialMap::coff(int i) {
	// TODO
	return m_Polynomial[i];
}

void PolynomialMap::compress() {
	// TODO
	PolynomialMap tmp(*this);
	m_Polynomial.clear();
	for (const auto& term : tmp.m_Polynomial) {
		if (fabs(term.second) > 1e-10) {
			coff(term.first) = term.second;
		}
	}
}

PolynomialMap PolynomialMap::operator+(const PolynomialMap& right) const {
	// TODO
	PolynomialMap tmp(right);
	for (const auto& term : m_Polynomial) {
		tmp.coff(term.first) += term.second;
	}
	tmp.compress();
	return tmp;
}

PolynomialMap PolynomialMap::operator-(const PolynomialMap& right) const {
	// TODO
	PolynomialMap tmp(right);
	for (const auto& term : m_Polynomial) {
		tmp.coff(term.first) -= term.second;
	}
	tmp.compress();
	return tmp;
}

PolynomialMap PolynomialMap::operator*(const PolynomialMap& right) const {
	// TODO
	PolynomialMap tmp = PolynomialMap();
	for (const auto& term : m_Polynomial) {
		for (const auto& term2 : right.m_Polynomial) {
			int deg = term.first + term2.first;
			double cof = term.second * term2.second;
			tmp.coff(deg) += cof;
		}
	}
	tmp.compress();
	return tmp;
}

PolynomialMap& PolynomialMap::operator=(const PolynomialMap& right) {
	// TODO
	m_Polynomial = right.m_Polynomial;
	return *this;
}

void PolynomialMap::Print() const {
	// TODO
	auto itr = m_Polynomial.begin();
	while (itr != m_Polynomial.end()) {
		if (itr != m_Polynomial.begin()) {
			if (itr->second > 0)
				cout << "+";
		}
		if (!abs(itr->second - 1) < 1e-10) {
			cout << itr->second;
		}
		if (itr->first > 0)
			cout << "x^" << itr->second;
		itr++;
	}
	cout << endl;
}

bool PolynomialMap::ReadFromFile(const string& file) {
	m_Polynomial.clear();
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
		m_Polynomial.emplace(deg, cof);
	}
	fin.close();
	return true;
}
