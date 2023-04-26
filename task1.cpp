#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>
using namespace std;
namespace MatrixUtils {
    class Matrix {
    public:
        int n;
        int m;
        vector<vector<double>> data;

        Matrix(int n, int m) {
            this->n = n;
            this->m = m;
            this->data = *new vector<vector<double>>(n, vector<double>(m));
        }

        Matrix operator+(Matrix b) const {
            if (this->n != b.n || this->m != b.m) {
                cout << "Error: the dimensional problem occurred" << endl;
                Matrix result = *new Matrix(0, 0);
                return result;
            }
            Matrix result = *new Matrix(n, m);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < m; j++) {
                    result.data[i].at(j) = this->data[i].at(j) + b.data[i].at(j);
                }
            }
            return result;
        }
        Matrix operator-(Matrix b) const {
            if (this->n != b.n || this->m != b.m) {
                cout << "Error: the dimensional problem occurred" << endl;
                Matrix result = *new Matrix(0, 0);
                return result;
            }
            Matrix result = *new Matrix(n, m);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < m; j++) {
                    result.data[i].at(j) = this->data[i].at(j) - b.data[i].at(j);
                }
            }
            return result;
        }
        Matrix operator*(Matrix b) const {
            if (this->m != b.n) {
                cout << "Error: the dimensional problem occurred" << endl;
                Matrix result = *new Matrix(0, 0);
                return result;
            }
            Matrix result = *new Matrix(this->n, b.m);
            for (int i = 0; i < result.n; i++) {
                for (int j = 0; j < result.m; j++) {
                    double sum = 0;
                    for (int k = 0; k < b.n; k++) {
                        sum += this->data[i].at(k) * b.data[k].at(j);
                    }
                    result.data[i].at(j) = sum;
                }
            }
            return result;
        }
        Matrix operator!() const {
            Matrix result = *new Matrix(this->m, this->n);
            for (int i = 0; i < result.n; i++) {
                for (int j = 0; j < result.m; j++) {
                    result.data[i].at(j) = this->data[j].at(i);
                }
            }
            return result;
        }
        void operator=(Matrix b) {
            this->n = b.n;
            this->m = b.m;
            this->data = b.data;
        }

        friend ostream& operator<<(ostream& out, Matrix mtrx) {
            for (int i = 0; i < mtrx.n; i++) {
                for (int j = 0; j < mtrx.m; j++) {
                    cout.precision(4);
                    cout << fixed << mtrx.data[i].at(j) << " ";
                }
                cout << endl;
            }
            return cout;
        };

        friend Matrix operator<<(Matrix mtrx1, Matrix mtrx2) {
            int M_new = mtrx1.m + mtrx2.m;
            Matrix composed = *new Matrix(mtrx1.n, M_new);
            for (int i = 0; i < mtrx1.n; i++) {
                int k = 0;
                for (int j = 0; j < mtrx1.m; j++) {
                    composed.data[i].at(j) = mtrx1.data[i].at(j);
                    k++;
                }
                for (int j = 0; j < mtrx2.m; j++) {
                    composed.data[i].at(k) = mtrx2.data[i].at(j);
                    k++;
                }
            }
            return composed;
        }

        friend istream& operator>>(istream& in, Matrix& mtrx) {
            string num;
            for (int i = 0; i < mtrx.n; i++) {
                for (int j = 0; j < mtrx.m; j++) {
                    cin >> num;
                    mtrx.data[i].at(j) = (stod(num));
                }
            }
            return cin;
        };
    };
    class SquareMatrix : public Matrix {
    public:
        SquareMatrix(int n);
    };
    SquareMatrix::SquareMatrix(int n) : Matrix(n, n) {}

    class IdentityMatrix : public SquareMatrix {
    public:
        IdentityMatrix(int n);
    };
    IdentityMatrix::IdentityMatrix(int n) : SquareMatrix(n) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                this->data[i].at(j) = ((i == j) ? 1 : 0);
            }
        }
    };

    class EliminationMatrix : public IdentityMatrix {
    public:
        EliminationMatrix(int n, int i, int j, Matrix A);
    };
    EliminationMatrix::EliminationMatrix(int n, int i, int j, MatrixUtils::Matrix A) : IdentityMatrix(n) {
        double ost = double(A.data[j].at(i)) / A.data[i].at(i);
        for (int k = 0; k < n; k++) {
            this->data[j].at(k) -= (this->data[i].at(k) * (ost));
        }
    }

    class PermutationMatrix : public IdentityMatrix {
    public:
        PermutationMatrix(int n, int i, int j);
    };
    PermutationMatrix::PermutationMatrix(int n, int i, int j) : IdentityMatrix(n) {
        for (int k = 0; k < n; k++) {
            double temp = this->data[i].at(k);
            this->data[i].at(k) = this->data[j].at(k);
            this->data[j].at(k) = temp;
        }
    };

    class AugmentedMatrix : public Matrix {
    public:
        AugmentedMatrix(Matrix m1, Matrix m2);
    };
    AugmentedMatrix::AugmentedMatrix(MatrixUtils::Matrix m1, MatrixUtils::Matrix m2) : Matrix(m1.n, (m1.m + m2.m)) {
        Matrix aug = m1 << m2;
        this->data = aug.data;
    }
}
namespace VectorUtils {
    static double* findMax(vector<vector<double>> vector, int n, int j, int i) {
        static double max[2];
        *max = 0, *(max+1) = 0;
        for (i; i < n; i++) {
            if (abs(vector[i].at(j)) > abs(*max) ) {
                *max = vector[i].at(j);
                *(max+1) = i;
            }
        }
        return max;
    }
}

using namespace MatrixUtils;
#ifdef WIN32
#define GNUPLOT_NAME "C:\\gnuplot\\bin\\gnuplot -persist"
#else
#define GNUPLOT_NAME "gnuplot -persist"
#endif
int main() {
#ifdef WIN32
    FILE* pipe = _popen(GNUPLOT_NAME, "w");
#else
    FILE* pipe = popen(GNUPLOT_NAME, "w");
#endif
    int m;
    cin >> m;

    Matrix x = *new Matrix(m, 1);
    Matrix y = *new Matrix(m, 1);
    for (int i = 0; i < m; i++) {
        double x_i, y_i;
        cin >> x_i;
        cin >> y_i;
        x.data[i].at(0) = x_i;
        y.data[i].at(0) = y_i;
    }
    int n;
    cin >> n;
    Matrix A = *new Matrix(m, n+1);
    for (int i = 0; i < m; i++) {
        for (int j = n; j >= 0; j--) {
            A.data[i].at(j) = pow(x.data[i].at(0), j);
        }
    }
    Matrix A_T = !A;
    cout << "A:" << endl;
    cout << A;

    Matrix A_T_A = A_T * A;
    cout << "A_T*A:" << endl;
    cout << A_T_A;

    Matrix I = *new IdentityMatrix(A_T_A.n);
    Matrix aug = *new AugmentedMatrix(A_T_A, I);
    // Direct
    for (int i = 0; i < A_T_A.n-1; i++) {
        for (int j = i; j < A_T_A.n; j++) {
            if (i == j) {
                double* max = VectorUtils::findMax(aug.data, aug.n, j, i);
                if (abs(*max) > abs(aug.data[i].at(j)) && *(max+1) != i) {
                    // permutation
                    Matrix P = *new PermutationMatrix(A_T_A.n, i, (int)*(max+1));
                    aug = P * aug;
                }
                continue;
            }
            if (j != 0 && aug.data[j].at(i) != 0) {
                // elimination
                Matrix E = *new EliminationMatrix(A_T_A.n, i, j, aug);
                aug = E * aug;
            }
        }
    }
    // Way back
    for (int i = A_T_A.n-1; i > 0; i--) {
        for (int j = i-1; j >= 0; j--) {
            if (abs(aug.data[j].at(i)) != 0) { // != 0
                // elimination
                Matrix E = *new EliminationMatrix(A_T_A.n, i, j, aug);
                aug = E * aug;
            }
        }
    }
    // Diagonal normalization
    for (int i = A_T_A.n-1; i >= 0; i--) {
        double piv = aug.data[i].at(i);
        for (int j = 0; j < aug.m; j++) {
            aug.data[i].at(j) /= piv;
        }
    }

    Matrix A_T_A_inv = *new SquareMatrix(A_T_A.n);
    for (int i = 0; i < A_T_A.n; i++) {
        for (int j = 0; j < A_T_A.n; j++) {
            A_T_A_inv.data[i].at(j) = aug.data[i].at(j+A_T_A.n);
        }
    }

    cout << "(A_T*A)^-1:" << endl;
    cout << A_T_A_inv;

    Matrix A_T_b = A_T * y;
    cout << "A_T*b:" << endl;
    cout << A_T_b;

    cout << "x~:" << endl;
    Matrix xapp = A_T_A_inv * A_T_b;
    cout << xapp;

    fprintf(pipe, "plot [-15 : 15] [-40 : 40] %lf*x**3 + %lf*x**2 + %lf*x**1 + %lf*x**0 , '-' using 1:2 with points\n", xapp.data[3].at(0), xapp.data[2].at(0), xapp.data[1].at(0), xapp.data[0].at(0));
    for (int i = 0; i < m; i++) {
        fprintf(pipe, "%f\t%f\n", x.data[i].at(0), y.data[i].at(0));
    }
    fprintf(pipe, "e\n");
    fflush(pipe);

#ifdef WIN32
    _pclose(pipe);
#else
    pclose(pipe);
#endif
}