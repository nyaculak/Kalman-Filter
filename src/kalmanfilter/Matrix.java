package kalmanfilter; 

import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.util.Locale;

public class Matrix implements Cloneable, java.io.Serializable {
  private static final long serialVersionUID = 1L;

  private final double[][] data;

  private final int numRows, numCols;

  public Matrix(int m, int n) {
    this.numRows = m;
    this.numCols = n;
    data = new double[m][n];
  }

  public Matrix(int m, int n, double s) {
    this.numRows = m;
    this.numCols = n;
    data = new double[m][n];
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        data[i][j] = s;
      }
    }
  }

  public Matrix(double[][] A) {
    numRows = A.length;
    numCols = A[0].length;
    for (int i = 0; i < numRows; i++) {
      if (A[i].length != numCols) {
        throw new IllegalArgumentException(
            "All rows must have the same length.");
      }
    }
    this.data = A;
  }

  /**
   * Construct a matrix from a copy of a 2-D array.
   * 
   * @param A
   *            Two-dimensional array of doubles.
   * @exception IllegalArgumentException
   *                All rows must have the same length
   */

  public static Matrix constructWithCopy(double[][] A) {
    int m = A.length;
    int n = A[0].length;
    Matrix X = new Matrix(m, n);
    double[][] C = X.getArray();
    for (int i = 0; i < m; i++) {
      if (A[i].length != n) {
        throw new IllegalArgumentException(
            "All rows must have the same length.");
      }
      for (int j = 0; j < n; j++) {
        C[i][j] = A[i][j];
      }
    }
    return X;
  }

  /**
   * Make a deep copy of a matrix
   */

  public Matrix copy() {
    Matrix X = new Matrix(numRows, numCols);
    double[][] C = X.getArray();
    for (int i = 0; i < numRows; i++) {
      for (int j = 0; j < numCols; j++) {
        C[i][j] = data[i][j];
      }
    }
    return X;
  }

  /**
   * Clone the Matrix object.
   */

  @Override
  public Object clone() {
    return this.copy();
  }

  public double[][] getArray() {
    return data;
  }

  /**
   * Copy the internal two-dimensional array.
   * 
   * @return Two-dimensional array copy of matrix elements.
   */

  public double[][] getArrayCopy() {
    double[][] C = new double[numRows][numCols];
    for (int i = 0; i < numRows; i++) {
      for (int j = 0; j < numCols; j++) {
        C[i][j] = data[i][j];
      }
    }
    return C;
  }

  /**
   * Get row dimension.
   * 
   * @return m, the number of rows.
   */

  public int getRowDimension() {
    return numRows;
  }

  /**
   * Get column dimension.
   * 
   * @return n, the number of columns.
   */

  public int getColumnDimension() {
    return numCols;
  }

  /**
   * Get a single element.
   * 
   * @param i
   *            Row index.
   * @param j
   *            Column index.
   * @return A(i,j)
   * @exception ArrayIndexOutOfBoundsException
   */

  public double get(int i, int j) {
    return data[i][j];
  }

  /**
   * Get a submatrix.
   * 
   * @param i0
   *            Initial row index
   * @param i1
   *            Final row index
   * @param j0
   *            Initial column index
   * @param j1
   *            Final column index
   * @return A(i0:i1,j0:j1)
   * @exception ArrayIndexOutOfBoundsException
   *                Submatrix indices
   */

  public Matrix getMatrix(int i0, int i1, int j0, int j1) {
    Matrix X = new Matrix(i1 - i0 + 1, j1 - j0 + 1);
    double[][] B = X.getArray();
    try {
      for (int i = i0; i <= i1; i++) {
        for (int j = j0; j <= j1; j++) {
          B[i - i0][j - j0] = data[i][j];
        }
      }
    } catch (ArrayIndexOutOfBoundsException e) {
      throw new ArrayIndexOutOfBoundsException("Submatrix indices");
    }
    return X;
  }


  /**
   * Get a submatrix.
   * 
   * @param r
   *            Array of row indices.
   * @param j0
   *            Initial column index
   * @param j1
   *            Final column index
   * @return A(r(:),j0:j1)
   * @exception ArrayIndexOutOfBoundsException
   *                Submatrix indices
   */

  public Matrix getMatrix(int[] r, int j0, int j1) {
    Matrix X = new Matrix(r.length, j1 - j0 + 1);
    double[][] B = X.getArray();
    try {
      for (int i = 0; i < r.length; i++) {
        for (int j = j0; j <= j1; j++) {
          B[i][j - j0] = data[r[i]][j];
        }
      }
    } catch (ArrayIndexOutOfBoundsException e) {
      throw new ArrayIndexOutOfBoundsException("Submatrix indices");
    }
    return X;
  }

  /**
   * Set a single element.
   * 
   * @param i
   *            Row index.
   * @param j
   *            Column index.
   * @param s
   *            A(i,j).
   * @exception ArrayIndexOutOfBoundsException
   */

  public void set(int i, int j, double s) {
    data[i][j] = s;
  }

  /**
   * Set a submatrix.
   * 
   * @param i0
   *            Initial row index
   * @param i1
   *            Final row index
   * @param j0
   *            Initial column index
   * @param j1
   *            Final column index
   * @param X
   *            A(i0:i1,j0:j1)
   * @exception ArrayIndexOutOfBoundsException
   *                Submatrix indices
   */

  public void setMatrix(int i0, int i1, int j0, int j1, Matrix X) {
    try {
      for (int i = i0; i <= i1; i++) {
        for (int j = j0; j <= j1; j++) {
          data[i][j] = X.get(i - i0, j - j0);
        }
      }
    } catch (ArrayIndexOutOfBoundsException e) {
      throw new ArrayIndexOutOfBoundsException("Submatrix indices");
    }
  }

  /**
   * Set a submatrix.
   * 
   * @param r
   *            Array of row indices.
   * @param c
   *            Array of column indices.
   * @param X
   *            A(r(:),c(:))
   * @exception ArrayIndexOutOfBoundsException
   *                Submatrix indices
   */

  public void setMatrix(int[] r, int[] c, Matrix X) {
    try {
      for (int i = 0; i < r.length; i++) {
        for (int j = 0; j < c.length; j++) {
          data[r[i]][c[j]] = X.get(i, j);
        }
      }
    } catch (ArrayIndexOutOfBoundsException e) {
      throw new ArrayIndexOutOfBoundsException("Submatrix indices");
    }
  }

  public Matrix transpose() {
    Matrix X = new Matrix(numCols, numRows);
    double[][] C = X.getArray();
    for (int i = 0; i < numRows; i++) {
      for (int j = 0; j < numCols; j++) {
        C[j][i] = data[i][j];
      }
    }
    return X;
  }

  public Matrix plus(Matrix B) {
    checkMatrixDimensions(B);
    Matrix X = new Matrix(numRows, numCols);
    double[][] C = X.getArray();
    for (int i = 0; i < numRows; i++) {
      for (int j = 0; j < numCols; j++) {
        C[i][j] = data[i][j] + B.data[i][j];
      }
    }
    return X;
  }

  public Matrix plusEquals(Matrix B) {
    checkMatrixDimensions(B);
    for (int i = 0; i < numRows; i++) {
      for (int j = 0; j < numCols; j++) {
        data[i][j] = data[i][j] + B.data[i][j];
      }
    }
    return this;
  }

  public Matrix minus(Matrix B) {
    checkMatrixDimensions(B);
    Matrix X = new Matrix(numRows, numCols);
    double[][] C = X.getArray();
    for (int i = 0; i < numRows; i++) {
      for (int j = 0; j < numCols; j++) {
        C[i][j] = data[i][j] - B.data[i][j];
      }
    }
    return X;
  }

  public Matrix minusEquals(Matrix B) {
    checkMatrixDimensions(B);
    for (int i = 0; i < numRows; i++) {
      for (int j = 0; j < numCols; j++) {
        data[i][j] = data[i][j] - B.data[i][j];
      }
    }
    return this;
  }

  /**
   * Multiply a matrix by a scalar, C = s*A
   * 
   * @param s
   *            scalar
   * @return s*A
   */

  public Matrix times(double s) {
    Matrix X = new Matrix(numRows, numCols);
    double[][] C = X.getArray();
    for (int i = 0; i < numRows; i++) {
      for (int j = 0; j < numCols; j++) {
        C[i][j] = s * data[i][j];
      }
    }
    return X;
  }

  /**
   * Multiply a matrix by a scalar in place, A = s*A
   * 
   * @param s scalar
   * @return replace A by s*A
   */

  public Matrix timesEquals(double s) {
    for (int i = 0; i < numRows; i++) {
      for (int j = 0; j < numCols; j++) {
        data[i][j] = s * data[i][j];
      }
    }
    return this;
  }

  /**
   * Linear algebraic matrix multiplication, A * B
   * 
   * @param B
   *            another matrix
   * @return Matrix product, A * B
   * @exception IllegalArgumentException
   *                Matrix inner dimensions must agree.
   */

  public Matrix times(Matrix B) {
    if (B.numRows != numCols) {
      throw new IllegalArgumentException(
          "Matrix inner dimensions must agree.");
    }
    Matrix X = new Matrix(numRows, B.numCols);
    double[][] C = X.getArray();
    double[] Bcolj = new double[numCols];
    for (int j = 0; j < B.numCols; j++) {
      for (int k = 0; k < numCols; k++) {
        Bcolj[k] = B.data[k][j];
      }
      for (int i = 0; i < numRows; i++) {
        double[] Arowi = data[i];
        double s = 0;
        for (int k = 0; k < numCols; k++) {
          s += Arowi[k] * Bcolj[k];
        }
        C[i][j] = s;
      }
    }
    return X;
  }

  /**
   * Solve A*X = B
   * 
   * @param B
   *            right hand side
   * @return solution if A is square, least squares solution otherwise
   */

  public Matrix solve(Matrix B) {
    // assumed m == n
    return new LUDecomposition(this).solve(B);

  }

  /**
   * Solve X*A = B, which is also A'*X' = B'
   * 
   * @param B
   *            right hand side
   * @return solution if A is square, least squares solution otherwise.
   */

  public Matrix solveTranspose(Matrix B) {
    return transpose().solve(B.transpose());
  }

  /**
   * Matrix inverse or pseudoinverse
   * 
   * @return inverse(A) if A is square, pseudoinverse otherwise.
   */

  public Matrix inverse() {
    return solve(identity(numRows, numRows));
  }

  /**
   * Matrix determinant
   * 
   * @return determinant
   */

  public double det() {
    return new LUDecomposition(this).det();
  }

  /**
   * Generate identity matrix
   * 
   * @param m
   *            Number of rows.
   * @param n
   *            Number of colums.
   * @return An m-by-n matrix with ones on the diagonal and zeros elsewhere.
   */

  public static Matrix identity(int m, int n) {
    Matrix A = new Matrix(m, n);
    double[][] X = A.getArray();
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        X[i][j] = (i == j ? 1.0 : 0.0);
      }
    }
    return A;
  }

  /**
   * Print the matrix to stdout. Line the elements up in columns with a
   * Fortran-like 'Fw.d' style format.
   * 
   * @param w
   *            Column width.
   * @param d
   *            Number of digits after the decimal.
   */

  public void print(int w, int d) {
    print(new PrintWriter(System.out, true), w, d);
  }

  /**
   * Print the matrix to the output stream. Line the elements up in columns
   * with a Fortran-like 'Fw.d' style format.
   * 
   * @param output
   *            Output stream.
   * @param w
   *            Column width.
   * @param d
   *            Number of digits after the decimal.
   */

  public void print(PrintWriter output, int w, int d) {
    DecimalFormat format = new DecimalFormat();
    format.setDecimalFormatSymbols(new DecimalFormatSymbols(Locale.US));
    format.setMinimumIntegerDigits(1);
    format.setMaximumFractionDigits(d);
    format.setMinimumFractionDigits(d);
    format.setGroupingUsed(false);
    print(output, format, w + 2);
  }

  /**
   * Print the matrix to the output stream. Line the elements up in columns.
   * Use the format object, and right justify within columns of width
   * characters. Note that is the matrix is to be read back in, you probably
   * will want to use a NumberFormat that is set to US Locale.
   * 
   * @param output
   *            the output stream.
   * @param format
   *            A formatting object to format the matrix elements
   * @param width
   *            Column width.
   * @see java.text.DecimalFormat#setDecimalFormatSymbols
   */

  public void print(PrintWriter output, NumberFormat format, int width) {
    output.println(); // start on new line.
    for (int i = 0; i < numRows; i++) {
      for (int j = 0; j < numCols; j++) {
        String s = format.format(data[i][j]); // format the number
        int padding = Math.max(1, width - s.length()); // At _least_ 1
        // space
        for (int k = 0; k < padding; k++)
          output.print(' ');
        output.print(s);
      }
      output.println();
    }
    output.println(); // end with blank line.
  }

  @Override
  public String toString() {
    StringBuffer buf = new StringBuffer();
    for (int i = 0; i < getRowDimension(); i++) {

      for (int j = 0; j < getColumnDimension(); j++) {
        buf.append(get(i, j));
        buf.append(" ");
      }
      buf.append("\n");
    }

    return buf.toString();
  }

  /*
   * ------------------------ Private Methods ------------------------
   */

  /** Check if size(A) == size(B) * */

  private void checkMatrixDimensions(Matrix B) {
    if (B.numRows != numRows || B.numCols != numCols) {
      throw new IllegalArgumentException("Matrix dimensions must agree.");
    }
  }
}


/**
 * LU Decomposition.
 * <P>
 * For an m-by-n matrix A with m >= n, the LU decomposition is an m-by-n unit
 * lower triangular matrix L, an n-by-n upper triangular matrix U, and a
 * permutation vector piv of length m so that A(piv,:) = L*U. If m < n, then L
 * is m-by-m and U is m-by-n.
 * <P>
 * The LU decompostion with pivoting always exists, even if the matrix is
 * singular, so the constructor will never fail. The primary use of the LU
 * decomposition is in the solution of square systems of simultaneous linear
 * equations. This will fail if isNonsingular() returns false.
 */
 class LUDecomposition implements java.io.Serializable {
  private static final long serialVersionUID = 1L;

  /*
   * ------------------------ Class variables ------------------------
   */

  /**
   * Array for internal storage of decomposition.
   * 
   * @serial internal array storage.
   */
  private final double[][] LU;

  /**
   * Row and column dimensions, and pivot sign.
   * 
   * @serial column dimension.
   * @serial row dimension.
   * @serial pivot sign.
   */
  private final int m, n;

  private int pivsign;

  /**
   * Internal storage of pivot vector.
   * 
   * @serial pivot vector.
   */
  private final int[] piv;

  /*
   * ------------------------ Constructor ------------------------
   */

  /**
   * LU Decomposition, a structure to access L, U and piv.
   * 
   * @param A
   *            Rectangular matrix
   */
  public LUDecomposition(Matrix A) {

    // Use a "left-looking", dot-product, Crout/Doolittle algorithm.

    LU = A.getArrayCopy();
    m = A.getRowDimension();
    n = A.getColumnDimension();
    piv = new int[m];
    for (int i = 0; i < m; i++) {
      piv[i] = i;
    }
    pivsign = 1;
    double[] LUrowi;
    double[] LUcolj = new double[m];

    // Outer loop.

    for (int j = 0; j < n; j++) {

      // Make a copy of the j-th column to localize references.

      for (int i = 0; i < m; i++) {
        LUcolj[i] = LU[i][j];
      }

      // Apply previous transformations.

      for (int i = 0; i < m; i++) {
        LUrowi = LU[i];

        // Most of the time is spent in the following dot product.

        int kmax = Math.min(i, j);
        double s = 0.0;
        for (int k = 0; k < kmax; k++) {
          s += LUrowi[k] * LUcolj[k];
        }

        LUrowi[j] = LUcolj[i] -= s;
      }

      // Find pivot and exchange if necessary.

      int p = j;
      for (int i = j + 1; i < m; i++) {
        if (Math.abs(LUcolj[i]) > Math.abs(LUcolj[p])) {
          p = i;
        }
      }
      if (p != j) {
        for (int k = 0; k < n; k++) {
          double t = LU[p][k];
          LU[p][k] = LU[j][k];
          LU[j][k] = t;
        }
        int k = piv[p];
        piv[p] = piv[j];
        piv[j] = k;
        pivsign = -pivsign;
      }

      // Compute multipliers.

      if (j < m & LU[j][j] != 0.0) {
        for (int i = j + 1; i < m; i++) {
          LU[i][j] /= LU[j][j];
        }
      }
    }
  }

  /**
   * Is the matrix nonsingular?
   * 
   * @return true if U, and hence A, is nonsingular.
   */
  public boolean isNonsingular() {
    for (int j = 0; j < n; j++) {
      if (LU[j][j] == 0)
        return false;
    }
    return true;
  }

  /**
   * Determinant
   * 
   * @return det(A)
   * @exception IllegalArgumentException
   *                Matrix must be square
   */
  public double det() {
    if (m != n) {
      throw new IllegalArgumentException("Matrix must be square.");
    }
    double d = pivsign;
    for (int j = 0; j < n; j++) {
      d *= LU[j][j];
    }
    return d;
  }

  /**
   * Solve A*X = B
   * 
   * @param B
   *            A Matrix with as many rows as A and any number of columns.
   * @return X so that L*U*X = B(piv,:)
   * @exception IllegalArgumentException
   *                Matrix row dimensions must agree.
   * @exception RuntimeException
   *                Matrix is singular.
   */
  public Matrix solve(Matrix B) {
    if (B.getRowDimension() != m) {
      throw new IllegalArgumentException(
          "Matrix row dimensions must agree.");
    }
    if (!this.isNonsingular()) {
      throw new RuntimeException("Matrix is singular.");
    }

    // Copy right hand side with pivoting
    int nx = B.getColumnDimension();
    Matrix Xmat = B.getMatrix(piv, 0, nx - 1);
    double[][] X = Xmat.getArray();

    // Solve L*Y = B(piv,:)
    for (int k = 0; k < n; k++) {
      for (int i = k + 1; i < n; i++) {
        for (int j = 0; j < nx; j++) {
          X[i][j] -= X[k][j] * LU[i][k];
        }
      }
    }
    // Solve U*X = Y;
    for (int k = n - 1; k >= 0; k--) {
      for (int j = 0; j < nx; j++) {
        X[k][j] /= LU[k][k];
      }
      for (int i = 0; i < k; i++) {
        for (int j = 0; j < nx; j++) {
          X[i][j] -= X[k][j] * LU[i][k];
        }
      }
    }
    return Xmat;
  }
}