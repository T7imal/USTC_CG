#pragma once
#include <iostream>
template <typename T>
// interfaces of Dynamic Array class DArray
class DArray {
public:

	// default constructor
	DArray() {
		Init();
	}

	// set an array with default values
	DArray(int nSize, T dValue) {
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

	DArray(const DArray& arr) {
		//TODO
		m_nSize = arr.m_nSize;
		m_nMax = arr.m_nMax;
		m_pData = new T[m_nMax];
		for (int i = 0; i < m_nSize; ++i) {
			m_pData[i] = arr[i];
		}
	}

	// deconstructor
	~DArray() {
		Free();
	}

	// display the elements of the array
	void Print() const {
		//TODO
		for (int i = 0; i < m_nSize; ++i) {
			std::cout << m_pData[i] << " ";
		}
		std::cout << std::endl;
	}

	// initilize the array
	void Init() {
		//TODO
		m_nSize = 0;
		m_nMax = 1;
		m_pData = new T[1];
	}

	// free the array
	void Free() {
		//TODO
		delete(m_pData);
	}

	// get the size of the array
	int GetSize() const {
		//TODO
		return m_nSize;
	}

	// set the size of the array
	void SetSize(int nSize) {
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
	const T& GetAt(int nIndex) const {
		//TODO
		return m_pData[nIndex];
	}

	// set the value of an element 
	void SetAt(int nIndex, T dValue) {
		//TODO
		m_pData[nIndex] = dValue;
	}

	// overload operator '[]'
	T& operator[](int nIndex) {
		// TODO
		return m_pData[nIndex];
	}

	// overload operator '[]'
	const T& operator[](int nIndex) const {
		//TODO
		return m_pData[nIndex];
	}

	// add a new element at the end of the array
	void PushBack(T dValue) {
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
	void DeleteAt(int nIndex) {
		//TODO
		if (nIndex < 0 || nIndex >= m_nSize) return;
		for (int i = nIndex; i < m_nSize - 1; ++i) {
			m_pData[i] = m_pData[i + 1];
		}
		m_nSize--;
	}

	// insert a new element at some index
	void InsertAt(int nIndex, T dValue) {
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
	DArray& operator = (const DArray& arr) {
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


private:
	T* m_pData; // the pointer to the array memory
	int m_nSize; // the size of the array
	int m_nMax;

private:
	void Init(); // initilize the array
	void Free(); // free the array
	void Reserve(int nSize); // allocate enough memory
};
