#include <algorithm>
#include <random>
#include <vector>

#include <math.h>

struct Bounds {
	Bounds() = default;

	const Bounds& operator=(const Bounds& T) {
		this->m = T.m;
		this->k = T.k;

		return *this;
	}

	unsigned int m = 0 /*rows*/, k = 0 /*columns*/;
};

struct Matrix {
	Matrix() = default;

	void push(const float x) {
		val.push_back(x);
	}

	const Matrix& operator=(const Matrix& T) {
		this->val = T.val;
		this->bounds = T.bounds;

		return *this;
	}

	void operator+=(const Matrix& T);
	void operator-=(const Matrix& T);

	void operator*=(const Matrix& T);
	void operator*=(const float T);

	const Matrix operator+(const Matrix& T);
	const Matrix operator-(const Matrix& T);

	const Matrix operator*(const Matrix& T);
	const Matrix operator*(const float T);

	const Matrix Invert();

	std::vector<float> val;
	Bounds bounds;
};

void MoveMatrix(const Matrix& I, Matrix& T, const int offset_x = 0, const int offset_y = 0, const int offset_dx = 0, const int offset_dy = 0) {
	for (int i = offset_y; i < (I.bounds.m - offset_dy); i++) {
		for (int j = offset_x; j < (I.bounds.k - offset_dx); j++) {
			T.push(I.val[j + (i * I.bounds.k)]);
		}
	}

	T.bounds.m = I.bounds.m - (offset_y + offset_dy);
	T.bounds.k = I.bounds.k - (offset_x + offset_dx);
}

void MoveMatrixT(const Matrix& I, Matrix& T, const int offset_x = 0, const int offset_y = 0, const int offset_dx = 0, const int offset_dy = 0) {
	for (int i = offset_x; i < (I.bounds.k - offset_dx); i++) {
		for (int j = offset_y; j < (I.bounds.m - offset_dy); j++) {
			T.push(I.val[i + (j * I.bounds.k)]);
		}
	}

	T.bounds.m = I.bounds.k - (offset_y + offset_dy);
	T.bounds.k = I.bounds.m - (offset_x + offset_dx);
}

void ScalarMatrix(Matrix& I, const float scalar) {
	for (int i = 0; i < (I.bounds.m * I.bounds.k); i++)
		I.val[i] *= scalar;
}

void MulMatrix(Matrix& I, const Matrix& T) {
	if (I.bounds.k == T.bounds.k) {
		if (I.bounds.m == T.bounds.m) {

			for (int i = 0; i < (T.bounds.m * T.bounds.k); i++)
				I.val[i] *= T.val[i];

		}else if (T.bounds.m == 1) {

			for (int i = 0; i < I.bounds.m; i++)
				for (int j = 0; j < I.bounds.k; j++)
					I.val[j + (i * I.bounds.k)] *= T.val[j];

		}else if (I.bounds.m == 1) {
			Matrix T1;
			T1.bounds = T.bounds;

			for (int i = 0; i < T.bounds.m; i++)
				for (int j = 0; j < T.bounds.k; j++)
					T1.push(I.val[j] * T.val[j + (i * T.bounds.k)]);

			I = T1;
		}
	}else if (I.bounds.m == T.bounds.m) {
		if (T.bounds.k == 1) {

			for (int i = 0; i < I.bounds.k; i++)
				for (int j = 0; j < I.bounds.m; j++)
					I.val[i + (j * I.bounds.k)] *= T.val[j];

		}else if (I.bounds.k == 1) {
			Matrix T1;
			T1.bounds = T.bounds;

			for (int i = 0; i < T.bounds.m; i++)
				for (int j = 0; j < T.bounds.k; j++) {
					T1.push(I.val[i] * T.val[j + (i * T.bounds.k)]);
				}

			I = T1;
		}
	}
}

void AddMatrix(Matrix& I, const Matrix& T) {
	if (I.bounds.m == T.bounds.m) {
		if (I.bounds.k == T.bounds.k) {

			for (int i = 0; i < (T.bounds.k * T.bounds.m); i++)
				I.val[i] += T.val[i];

		}else if (T.bounds.k == 1) {

			for (int i = 0; i < I.bounds.m; i++)
				for (int j = 0; j < I.bounds.k; j++)
					I.val[j + (i * I.bounds.k)] += T.val[i];

		}else if (I.bounds.k == 1) {
			Matrix T1;
			T1.bounds = T.bounds;

			for (int i = 0; i < T.bounds.m; i++)
				for (int j = 0; j < T.bounds.k; j++)
					T1.push(I.val[i] + T.val[j + (i * T.bounds.k)]);

			I = T1;
		}
	}else if (I.bounds.k == T.bounds.k) {
		if (T.bounds.m == 1)

			for (int i = 0; i < I.bounds.m; i++)
				for (int j = 0; j < I.bounds.k; j++)
					I.val[j + (i * I.bounds.k)] += T.val[j];

		else if (I.bounds.m == 1) {
			Matrix T1;
			T1.bounds = T.bounds;

			for (int i = 0; i < T.bounds.m; i++)
				for (int j = 0; j < T.bounds.k; j++)
					T1.push(I.val[j] + T.val[j + (i * T.bounds.k)]);

			I = T1;
		}
	}
}

void DotMatrix(const Matrix& I, const Matrix& T, Matrix& R) {
	if (I.bounds.k == T.bounds.m) {
		if (T.bounds.k == 1) {
			for (int i = 0; i < I.bounds.m; i++) {
				float value = 0;

				for (int j = 0; j < I.bounds.k; j++) {
					value += I.val[j + (i * I.bounds.k)] * T.val[j];
				}

				R.push(value);
			}

			R.bounds.m = T.bounds.k;
			R.bounds.k = I.bounds.m;
		}else {
			for (int i = 0; i < I.bounds.m; i++) {
				for (int k = 0; k < T.bounds.k; k++) {
					float sum = 0;

					for (int j = 0; j < I.bounds.k; j++) {
						sum += I.val[j + (i * I.bounds.k)] * T.val[k + (j * T.bounds.k)];
					}

					R.push(sum);
				}
			}

			R.bounds.m = I.bounds.m;
			R.bounds.k = T.bounds.k;
		}
	}
}

void Matrix::operator+=(const Matrix& T) {
	AddMatrix(*this, T);
}

void Matrix::operator-=(const Matrix& T) {
	Matrix T1;
	T1 = T;
	T1.Invert();
	AddMatrix(*this, T1);
}

void Matrix::operator*=(const Matrix& T) {
	MulMatrix(*this, T);
}

void Matrix::operator*=(const float T) {
	ScalarMatrix(*this, T);
}

//

const Matrix Matrix::operator+(const Matrix& T) {
	Matrix C;
	C = *this;
	AddMatrix(C, T);

	return C;
}

const Matrix Matrix::operator-(const Matrix& T) {
	Matrix C, T1;
	C = *this;
	T1 = T;
	T1.Invert();
	AddMatrix(C, T1);

	return C;
}

const Matrix Matrix::operator*(const Matrix& T) {
	Matrix C;
	C = *this;
	MulMatrix(C, T);

	return C;
}

const Matrix Matrix::operator*(const float T) {
	Matrix C;
	C = *this;
	ScalarMatrix(C, T);

	return C;
}

const Matrix Matrix::Invert() {
	ScalarMatrix(*this, -1.0);
	return *this;
}