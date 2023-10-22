#pragma once
#include <iostream>
// interfaces of Dynamic Array class DArray
template <class T>
class DArray {
public:
	DArray(); // default constructor
	DArray(int nSize, T dValue = 0); // set an array with default values
	DArray(const DArray& arr); // copy constructor
	~DArray(); // deconstructor

	void Print() const; // print the elements of the array

	int GetSize() const; // get the size of the array
	void SetSize(int nSize); // set the size of the array

	const T& GetAt(int nIndex) const; // get an element at an index
	void SetAt(int nIndex, T dValue); // set the value of an element

	T& operator[](int nIndex); // overload operator '[]'
	const T& operator[](int nIndex) const; // overload operator '[]'

	void PushBack(T dValue); // add a new element at the end of the array
	void DeleteAt(int nIndex); // delete an element at some index
	void InsertAt(int nIndex, T dValue); // insert a new element at some index

	DArray& operator = (const DArray& arr); //overload operator '='

private:
	T* m_pData; // the pointer to the array memory
	int m_nSize; // the size of the array
	int m_nMax;

private:
	void Init(); // initilize the array
	void Free(); // free the array
	void Reserve(int nSize); // allocate enough memory
};

// default constructor
template<class T>
DArray<T>::DArray() {
	Init();
}

// set an array with default values
template<class T>
DArray<T>::DArray(int nSize, T dValue) {
	//TODO
	m_nSize = nSize;
	m_nMax = 1;
	while (m_nMax < nSize)
		m_nMax <<= 1;
	m_pData = new T[m_nMax];
	for (int i = 0; i < nSize; ++i) {
		m_pData[i] = dValue;
	}
}

template<class T>
DArray<T>::DArray(const DArray& arr) {
	//TODO
	m_nSize = arr.m_nSize;
	m_nMax = arr.m_nMax;
	m_pData = new T[m_nMax];
	for (int i = 0; i < m_nSize; ++i) {
		m_pData[i] = arr[i];
	}
}

// deconstructor
template<class T>
DArray<T>::~DArray() {
	Free();
}

// display the elements of the array
template<class T>
void DArray<T>::Print() const {
	//TODO
	for (int i = 0; i < m_nSize; ++i) {
		std::cout << m_pData[i] << " ";
	}
	std::cout << std::endl;
}

// initilize the array
template<class T>
void DArray<T>::Init() {
	//TODO
	m_nSize = 0;
	m_nMax = 1;
	m_pData = new T[1];
}

// free the array
template<class T>
void DArray<T>::Free() {
	//TODO
	delete(m_pData);
}

// get the size of the array
template<class T>
int DArray<T>::GetSize() const {
	//TODO
	return m_nSize;
}

// set the size of the array
template<class T>
void DArray<T>::SetSize(int nSize) {
	//TODO
	if (nSize == m_nSize || nSize < 0) return;
	if (nSize < m_nSize) {
		m_nSize = nSize;
	}
	else {
		if (nSize < m_nMax) {
			m_nSize = nSize;
			for (int i = m_nSize; i < nSize; ++i)
				m_pData[i] = 0;
		}
		else {
			while (m_nMax < nSize)
				m_nMax <<= 1;
			T* tmp = new T[m_nMax];
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
}

// get an element at an index
template<class T>
const T& DArray<T>::GetAt(int nIndex) const {
	//TODO
	return m_pData[nIndex];
}

// set the value of an element 
template<class T>
void DArray<T>::SetAt(int nIndex, T dValue) {
	//TODO
	m_pData[nIndex] = dValue;
}

// overload operator '[]'
template<class T>
T& DArray<T>::operator[](int nIndex) {
	// TODO
	return m_pData[nIndex];
}

// overload operator '[]'
template<class T>
const T& DArray<T>::operator[](int nIndex) const {
	//TODO
	return m_pData[nIndex];
}

// add a new element at the end of the array
template<class T>
void DArray<T>::PushBack(T dValue) {
	//TODO
	if (m_nSize < m_nMax)
		m_pData[m_nSize++] = dValue;
	else {
		m_nMax <<= 1;
		T* tmp = new T[m_nMax];
		for (int i = 0; i < m_nSize; ++i) {
			tmp[i] = m_pData[i];
		}
		tmp[m_nSize++] = dValue;
		delete(m_pData);
		m_pData = tmp;
	}
}

// delete an element at some index
template<class T>
void DArray<T>::DeleteAt(int nIndex) {
	//TODO
	if (nIndex < 0 || nIndex >= m_nSize) return;
	for (int i = nIndex; i < m_nSize - 1; ++i) {
		m_pData[i] = m_pData[i + 1];
	}
	m_nSize--;
}

// insert a new element at some index
template<class T>
void DArray<T>::InsertAt(int nIndex, T dValue) {
	//TODO
	if (nIndex < 0 || nIndex > m_nSize) return;
	if (m_nSize < m_nMax) {
		for (int i = m_nSize; i > nIndex; --i)
			m_pData[i] = m_pData[i - 1];
		m_pData[nIndex] = dValue;
		m_nSize++;
	}
	else {
		m_nMax <<= 1;
		T* tmp = new T[m_nMax];
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
}

// overload operator '='
template<class T>
DArray<T>& DArray<T>::operator = (const DArray& arr) {
	//TODO
	if (this == &arr) return *this;
	else {
		delete(m_pData);
		m_nSize = arr.m_nSize;
		m_pData = new T[m_nSize];
		for (int i = 0; i < m_nSize; ++i) {
			m_pData[i] = arr[i];
		}
		return *this;
	}
}