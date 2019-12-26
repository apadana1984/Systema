#pragma once
#include<iostream>
#include<vector>

// MIT License
// Author: Iman Khademi, December 2019
// http://imankhademi.com

using namespace std;

typedef vector<float> RowVector;
typedef vector<float> ColumnVector;
typedef vector<RowVector> Vector2D;

class Matrix
{
private:
	int rows = 1;
	int columns = 1;
	Vector2D data;
public:
	// constructors
	Matrix(int r=1, int c=1);
	Matrix(Matrix& m);
	Matrix(Vector2D a);
	~Matrix() {  };
	// accessors
	int getSize(int dim) { return ((dim == 0) ? (rows) : (columns)); };
	float& element(int r, int c);
	Vector2D& getData() { return data; };
	Matrix Transpose(void);
	Matrix Slice(int rb, int cb, int re=-1, int ce=-1);
	int findNonZeroRowElement(int i);
	bool isSquare() { return (rows == columns); };
	// elementary row operations
	void rSwap(int r1,int r2); // swap rows r1 and r2
	void rMult(float a, int r);// multiply row r1 by a
	void rAdd(int r1, int r2, float a); // add row 2 to row 1 => r1 := r1+r2
	// matrix operations
	friend Matrix operator+(Matrix m1,Matrix m2);
	friend Matrix operator-(Matrix m1, Matrix m2);
	friend Matrix operator*(float a, Matrix m);
	friend Matrix operator*(Matrix m1, Matrix m2);
	friend ColumnVector operator*(Matrix m, ColumnVector v);
	// display
	friend ostream& operator<<(ostream& out, Matrix& m);
};

class diagonalMatrix :public Matrix
{
	diagonalMatrix(RowVector vec);
};