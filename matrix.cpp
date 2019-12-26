#include<iostream>
#include<fstream>
#include<vector>

#include"matrix.h"

// MIT License
// Author: Iman Khademi, December 2019
// http://imankhademi.com

#define MIN(a,b) (a>b?b:a)
#define MAX(a,b) (a>b?a:b)

using namespace std;

enum matrixError { INVALID_ROW_COUNT, INVALID_COLUMN_COUNT , DIMENSIONS_NOT_MATCH};

Matrix diagonalMatrix(RowVector);
Matrix matAugment(Matrix, Matrix, int);
Matrix identityMatrix(int );
Matrix elementaryMatrix(Matrix );

Matrix::Matrix(int r, int c)
{
	try
	{
		if (r <= 0) throw(INVALID_ROW_COUNT);
		if (c <= 0) throw(INVALID_COLUMN_COUNT);

	}
	catch (matrixError error)
	{
		ofstream logFile("user-log.txt",ios::app);
		switch (error)
		{
			case INVALID_ROW_COUNT:
			{
				logFile << "Invalid Number of Rows:" << ", rows=" << r << endl;
				r = 1;
				break;
			}
			case INVALID_COLUMN_COUNT:
			{
				logFile << "Invalid Number of Columns:" << ", columns=" << c << endl;
				c = 1;
				break;
			}
		}
		logFile.close();
	}
	rows = r;
	columns = c;
	data = Vector2D(r, RowVector(c));
}

Matrix::Matrix(Matrix & m)
{
	rows = m.getSize(0);
	columns = m.getSize(1);
	data = Vector2D(rows,RowVector(columns));
	for (int i=0;i<rows;i++)
		for(int j=0;j<columns;j++)
		{
			data[i][j] = m.element(i,j);
		}
}

Matrix::Matrix(Vector2D a)
{
	rows = a.size();
	columns = a[0].size();
	data = Vector2D(rows, RowVector(columns));
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < columns; j++)
			data[i][j] = a[i][j];
}

Matrix Matrix::Transpose()
{
	Matrix mt(columns,rows);
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < columns; j++)
			mt.element(j,i) = data[i][j];
	return mt;
}

/* ********** Elementary Row Operations *************/

void Matrix::rSwap(int r1, int r2)
{
	RowVector vaux = data[r1];
	data[r1] = data[r2];
	data[r2] = vaux;
	return;
}

void Matrix::rMult(float a, int r)
{
	for (int j=0;j<columns;j++)	data[r][j] = a* data[r][j];
	return;
}
//r1 = r1+a*r2
void Matrix::rAdd(int r1, int r2, float a)
{
	for (int i = 0; i < columns;i++) data[r1][i] += a*data[r2][i];
	return;
}
/* ********************** Find the first zero element *****************************/
int findNonZero(RowVector vec)
{
	for (int i = 0; i < (int)vec.size(); i++)
		if (vec[i] != 0) return i;
	return -1;
}

/* *********************** Inverse and Elementary Form of a Matrix **********************/
Matrix Inverse(Matrix m)
{
	int rows = m.getSize(0);
	int columns = m.getSize(1);
	try
	{
		if (rows != columns) throw 1;
	}
	catch(...)
	{
		ofstream logFile("user-log.txt", ios::app);
		logFile << "Not a square matrix! Try 'pInverse'" << endl;
		logFile.close();
		return m;
	}
	Matrix I(rows,columns),MA(rows,2*columns);
	I = identityMatrix(m.getSize(0));
	MA = matAugment(m, I, 0); // [M | I ]
	MA = elementaryMatrix(MA);
	Matrix minv(rows,columns);
	minv = MA.Slice(0, m.getSize(0));
	return minv;
}

Matrix elementaryMatrix(Matrix m)
{
	Matrix me(m);
	int rows = me.getSize(0);
	int columns = me.getSize(1);
	int rowInd = 0;
	for (int i = 0; i < columns; i++) // column by column
	{
		int k = me.findNonZeroRowElement(i); // index of the row with non-zeros element in column 'i'
		if (k != -1)
		{
			if (k > rowInd) me.rSwap(k, rowInd);
			RowVector& curRow = me.getData()[rowInd];
			me.rMult(1.0f / curRow[i], rowInd);// set the first nonzero element to one
			for (int t = 0; t < rows; t++)
			{
				if (t != rowInd)
				{
					me.rAdd(t, rowInd, -me.getData()[t][i]);// set other elements on the column to zero
				}
			}
		}
		rowInd++;
	}
	return me;
}

int Matrix::findNonZeroRowElement(int col)// row index of the first non-zero element in a column
{
	for (int i = col; i < rows; i++)
		if (getData()[i][col] != 0) return i;
	return -1;
}


Matrix Matrix::Slice(int rb, int cb, int re, int ce)
{
	int ree=re, cee=ce;
	if (re == -1) ree = rows - 1;
	if (ce == -1) cee = columns - 1;
	try
	{
		if (rb >= rows || cb >= columns || re >= rows || ce >= columns) throw 1;
		if (rb > ree || cb > cee) throw 2;
	}
	catch (int err)
	{
		ofstream logFile("user-log.txt", ios::app);
		switch (err)
		{
		case 1:
			logFile << "Slice: Index out of bounds" << endl;
			break;
		case 2:
			logFile << "Slice: Upper index should be greater than the lower index" << endl;
			break;
		}
		logFile.close();
		Matrix m(rows, columns);
		return m;
	}
	Matrix m(ree - rb + 1, cee - cb + 1);
	for (int i = rb; i <= ree; i++)
		for (int j = cb; j <= cee; j++)
			m.element(i - rb, j - cb) = data[i][j];
	return m;
}

float & Matrix::element(int r, int c)
{
	if (r<rows && r >= 0 && c >= 0 && c<columns) return data[r][c];
	else return data[0][0];
}

vector<float> operator+(vector<float> v1, vector<float> v2)
{
	vector<float> v3;
	for (size_t i = 0; i < v1.size(); i++) v3.push_back(v1[i] + v2[i]);
	return v3;
}

vector<float> operator-(vector<float> v1, vector<float> v2)
{
	vector<float> v3;
	for (size_t i = 0; i < v1.size(); i++) v3.push_back(v1[i] - v2[i]);
	return v3;
}

Matrix operator+(Matrix m1, Matrix m2)
{
	int r1 = m1.getSize(0);
	int c1 = m1.getSize(1);
	try
	{
		if (r1 != m2.getSize(0) || c1 != m2.getSize(1)) throw(DIMENSIONS_NOT_MATCH);
	}
	catch(...)
	{
		ofstream logFile("user-log.txt", ios::app);
		logFile << "Dimensions don't match for matrix summation!" << endl;
		logFile.close();
		r1 = MIN(r1, m2.getSize(0));
		c1 = MIN(c1, m2.getSize(1));
	}
	Matrix m(m1);
	for (int i = 0; i < r1; i++)
		for (int j = 0; j < c1; j++)
		{
			m.element(i,j) += m2.element(i, j);
		}
	return m;
}

Matrix operator-(Matrix m1, Matrix m2)
{
	int r1 = m1.getSize(0);
	int c1 = m1.getSize(1);
	try
	{
		if (r1 != m2.getSize(0) || c1 != m2.getSize(1)) throw(DIMENSIONS_NOT_MATCH);
	}
	catch (...)
	{
		ofstream logFile("user-log.txt", ios::app);
		logFile << "Dimensions don't match for matrix summation!" << endl;
		logFile.close();
		r1 = MIN(r1, m2.getSize(0));
		c1 = MIN(c1, m2.getSize(1));
	}
	Matrix m(m1);
	for (int i = 0; i < r1; i++)
		for (int j = 0; j < c1; j++)
		{
			m.element(i, j) -= m2.element(i, j);
		}
	return m;
}
// m = m1*m2
Matrix operator*(Matrix m1, Matrix m2)
{
	int r2 = m2.getSize(0);
	int c1 = m1.getSize(1);
	try
	{
		if (r2 != c1) throw(DIMENSIONS_NOT_MATCH);
	}
	catch (...)
	{
		ofstream logFile("user-log.txt", ios::app);
		logFile << "Dimensions don't match for matrix multiplication!" << endl;
		logFile.close();
		return 0;
	}
	Matrix m(m1.getSize(0),m2.getSize(1));
	for (int i = 0; i < m1.getSize(0); i++)
		for (int j = 0; j < m2.getSize(1); j++)
		{
			m.element(i, j) = 0.f;
			for (int k = 0; k < r2; k++) m.element(i, j) += m1.element(i, k) * m2.element(k, j);
		}
	return m;
}

float innerProduct(vector<float>v1, vector<float>v2)
{
	float x = 0.f;
	for (size_t i = 0; i < v1.size(); i++) x += v1[i] * v2[i];
	return x;
}

// matrix-vector multiplication
ColumnVector operator*(Matrix m, ColumnVector v)
{
	try
	{
		if (m.getSize(1) != (int)v.size()) throw 1;
	}
	catch (...)
	{
		ofstream logFile("user-log.txt", ios::app);
		logFile << "Dimensions don't match for matrix-vector multiplication!" << endl;
		logFile.close();
		return v;
	}
	ColumnVector vo;
	for (RowVector vm : m.getData()) vo.push_back(innerProduct(vm,v));
	return vo;
}

// scalar-matrix multiplication: m1 = a*m

Matrix operator*(float a, Matrix m)
{
	Matrix m1(m);
	for (int i = 0; i < m1.getSize(0); i++)
		for (int j = 0; j < m1.getSize(1); j++)
		{
			m1.element(i, j) *= a;
		}
	return m1;
}

// scalar-vector multiplication
RowVector operator*(float a, RowVector v)
{
	RowVector va = v;
	for (float & x : va) x *= a;
	return va;
}

ostream & operator<<(ostream & out, Matrix & m)
{
	for (auto vec : m.data)
	{
		for (auto x : vec) out << x << ",";
		out << endl;
	}
	return out;
}

ostream& operator<<(ostream& out, RowVector v)
{
	for (float x : v)
	{
		out << x << ",";
	}
	return out;
}

ostream& operator<<(ostream& out, Vector2D v)
{
	for (RowVector vr : v)
	{
		out << vr << endl;
	}
	return out;
}
/* *************** Other Functions ************ */
Matrix matAugment(Matrix m1, Matrix m2, int aug = 0)
{
	try
	{
		if (m1.getSize(aug) != m2.getSize(aug)) throw 1;
	}
	catch (...)
	{
		ofstream logFile("user-log.txt", ios::app);
		logFile << "matAugment: Dimensions mismatch!";
		logFile.close();
		return m1;
	}
	Vector2D m1_data = m1.getData();
	Vector2D m2_data = m2.getData();
	if (aug == 0)
	{
		Vector2D ma_data = m1.getData();
		for (int i = 0 ; i < m1.getSize(0) ; i++ ) ma_data[i].insert(end(ma_data[i]), begin(m2_data[i]), end(m2_data[i]));
		Matrix ma(ma_data);
		return ma;
	}
	else
	{
		Vector2D ma_data = m1.getData();
		ma_data.insert(end(ma_data), begin(m2_data), end(m2_data));
		Matrix ma(ma_data);
		return ma;
	}
}

Vector2D Transpose(Vector2D A)
{
	size_t rows = A.size();
	size_t columns = A[0].size();
	Vector2D At(columns,RowVector(rows));
	for (size_t i = 0; i < rows; i++)
		for (size_t j = 0; j < columns; j++)
			At[j][i] = A[i][j];
	return At;
}

/************************ Special Matrices ************************/
Matrix diagonalMatrix(RowVector vec)
{
	size_t size = vec.size();
	Matrix m(size,size);
	for (size_t i = 0; i < size; i++) m.getData()[i][i] = vec[i];
	return m;
}

Matrix identityMatrix(int size)
{
	RowVector vec;
	vec.assign(size, 1.0f);
	Matrix m(size,size);
	m = diagonalMatrix(vec);
	return m;
}

/* **************** Matrix Decompositions ******************/
Vector2D operator+(Vector2D v1, Vector2D v2)
{
	Vector2D v(v1.size(),RowVector(v1[0].size()));
	for (size_t i = 0; i < v1.size(); i++)
		for (size_t j = 0; j < v1[0].size(); j++)
			v[i][j] = v1[i][j] + v2[i][j];
	return v;
}

Vector2D operator-(Vector2D v1, Vector2D v2)
{
	Vector2D v(v1.size(), RowVector(v1[0].size()));
	for (size_t i = 0; i < v1.size(); i++)
		for (size_t j = 0; j < v1[0].size(); j++)
			v[i][j] = v1[i][j] - v2[i][j];
	return v;
}

// scalar-vector2d multiplication
Vector2D operator*(float a, Vector2D v)
{
	Vector2D va = v;
	for (RowVector& vec : va) vec = a * vec;
	return va;
}

Vector2D operator*(Vector2D vl, Vector2D vr)
{
	Vector2D v(vl.size(),RowVector(vr[0].size()));
	for (size_t i = 0; i < vl.size(); i++)
		for (size_t j = 0; j < vr[0].size(); j++)
			for (size_t k = 0; k < vl[0].size(); k++)
				v[i][j] += vl[i][k] * vr[k][j];
	return v;
}

ColumnVector operator*(Vector2D vl, ColumnVector vr)
{
	ColumnVector v(vl.size());
	for (size_t i = 0; i < vl.size(); i++)
			for (size_t k = 0; k < vl[0].size(); k++)
				v[i] += vl[i][k] * vr[k];
	return v;
}


Vector2D Slice(Vector2D v, int rowb, int rowe, int colb, int cole)
{
	Vector2D vs(rowe - rowb + 1,RowVector(cole-colb+1));
	for (size_t i = rowb; i <= rowe; i++)
		for (size_t j = colb; j <= cole; j++)
			vs[i-rowb][j-colb] = v[i][j];
	return vs;
}

Vector2D getA21(RowVector v)
{
	Vector2D A21;
	for (size_t i = 1; i < v.size(); i++)
	{
		RowVector v1;
		v1.push_back(v[i]);
		A21.push_back(v1);
	}
	return A21;
}

Vector2D getA22(Vector2D A)
{
	Vector2D A22 = Slice(A, 1, A.size()-1, 1, A[0].size()-1);
	return A22;
}

void insertL21(Vector2D & L, Vector2D L21, int step)
{
	for (int i = 0; i < L21.size(); i++)
		L[step + i +1][step] = L21[i][0];
	return;
}

// Danger: Not checking weather M is positive semi-definite or not.
// Cholesky Factorization
// M = L*L' where M is positive semi-definite, L is lower triangular with positive diagonal elements
Matrix Cholesky(Matrix M)
{
	int n = M.getSize(0);
	Vector2D A = M.getData();
	Vector2D Lv(n, RowVector(n));
	for (int i = 0; i < n; i++)
	{
		Vector2D L21(n-i-1,RowVector(1));
		float l11 = sqrt(A[0][0]);
		if (A.size() > 1)
		{
			if (l11!=0) L21 = 1.0f / l11 * getA21(A[0]);
			A = getA22(A) - L21 * Transpose(L21);
			insertL21(Lv, L21, i);
		}
		Lv[i][i] = l11;
	}
	Matrix L(Lv);
	return Lv;
}