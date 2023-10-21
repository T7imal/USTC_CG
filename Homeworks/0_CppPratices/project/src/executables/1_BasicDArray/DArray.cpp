// implementation of class DArray
#include "DArray.h"
#include <iostream>

// default constructor
DArray::DArray() {
	Init();
}

// set an array with default values
DArray::DArray(int nSize, double dValue) {
	//TODO
	m_nSize = nSize;
	m_pData = new double[m_nSize];
	for (int i = 0; i < nSize; ++i) {
		m_pData[i] = dValue;
	}
}

DArray::DArray(const DArray& arr) {
	//TODO
	m_nSize = arr.m_nSize;
	m_pData = new double[m_nSize];
	for (int i = 0; i < m_nSize; ++i) {
		m_pData[i] = arr[i];
	}
}

// deconstructor
DArray::~DArray() {
	Free();
}

// display the elements of the array
void DArray::Print() const {
	//TODO
	for (int i = 0; i < m_nSize; ++i) {
		std::cout << m_pData[i] << " ";
	}
	std::cout << std::endl;
}

// initilize the array
void DArray::Init() {
	//TODO
	m_nSize = 0;
	m_pData = nullptr;
}

// free the array
void DArray::Free() {
	//TODO
	delete(m_pData);
}

// get the size of the array
int DArray::GetSize() const {
	//TODO
	return m_nSize;
}

// set the size of the array
void DArray::SetSize(int nSize) {
	//TODO
	if (nSize == m_nSize || nSize < 0) return;
	if (nSize < m_nSize) {
		double* tmp = new double[nSize];
		for (int i = 0; i < nSize; ++i) {
			tmp[i] = m_pData[i];
		}
		delete(m_pData);
		m_nSize = nSize;
		m_pData = tmp;
		return;
	}
	if (nSize > m_nSize) {
		double* tmp = new double[nSize];
		for (int i = 0; i < m_nSize; ++i) {
			tmp[i] = m_pData[i];
		}
		for (int i = m_nSize; i < nSize; ++i) {
			tmp[i] = 0;
		}
		delete(m_pData);
		m_nSize = nSize;
		m_pData = tmp;
	}
}

// get an element at an index
const double& DArray::GetAt(int nIndex) const {
	//TODO
	return m_pData[nIndex];
}

// set the value of an element 
void DArray::SetAt(int nIndex, double dValue) {
	//TODO
	m_pData[nIndex] = dValue;
}

// overload operator '[]'
const double& DArray::operator[](int nIndex) const {
	//TODO
	return m_pData[nIndex];
}

// add a new element at the end of the array
void DArray::PushBack(double dValue) {
	//TODO
	double* tmp = new double[m_nSize + 1];
	for (int i = 0; i < m_nSize; ++i) {
		tmp[i] = m_pData[i];
	}
	tmp[m_nSize] = dValue;
	delete(m_pData);
	m_nSize++;
	m_pData = tmp;
}

// delete an element at some index
void DArray::DeleteAt(int nIndex) {
	//TODO
	if(nIndex<0||nIndex>=m_nSize) return;
	double* tmp = new double[m_nSize - 1];
	for (int i = 0; i < nIndex; ++i) {
		tmp[i] = m_pData[i];
	}
	for (int i = nIndex; i < m_nSize - 1; ++i) {
		tmp[i] = m_pData[i + 1];
	}
	delete(m_pData);
	m_nSize--;
	m_pData = tmp;
}

// insert a new element at some index
void DArray::InsertAt(int nIndex, double dValue) {
	//TODO
	if (nIndex < 0 || nIndex > m_nSize) return;
	double* tmp = new double[m_nSize + 1];
	for (int i = 0; i < nIndex; ++i) {
		tmp[i] = m_pData[i];
	}
	tmp[nIndex] = dValue;
	for (int i = nIndex + 1; i < m_nSize + 1; ++i) {
		tmp[i] = m_pData[i - 1];
	}
	delete(m_pData);
	m_nSize++;
	m_pData = tmp;
}

// overload operator '='
DArray& DArray::operator = (const DArray& arr) {
	//TODO
	if (this == &arr) return *this;
	else {
		delete(m_pData);
		m_nSize = arr.m_nSize;
		m_pData = new double[m_nSize];
		for (int i = 0; i < m_nSize; ++i) {
			m_pData[i] = arr[i];
		}
		return *this;
	}
}
