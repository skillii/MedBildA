#ifndef NUMERICSHELPER_H
#define NUMERICSHELPER_H

#include "Definitions.h"
#include <cmath>

class NumericsHelper
{
public:

  static float min(const float a[], int length, int* min_index)
  {
      float min;
      int ind;

      min = a[0];
      ind = 0;

      for(int i = 1; i < length; i++)
      {
          if(a[i] < min)
          {
              min = a[0];
              ind = i;
          }
      }

      *min_index = ind;

      return min;
  }

  static float trilinearInterp(FloatImageType::Pointer inputImage, const float x_transformed, const float y_transformed, const float z_transformed)
  {
    // here starts the trilinear interpolation code

    // these are the interpolation weights
    const float xd = x_transformed - floor(x_transformed);
    const float yd = y_transformed - floor(y_transformed);
    const float zd = z_transformed - floor(z_transformed);

    const unsigned int input_depth = inputImage->GetLargestPossibleRegion().GetSize()[2];
    const unsigned int input_height = inputImage->GetLargestPossibleRegion().GetSize()[1];
    const unsigned int input_width = inputImage->GetLargestPossibleRegion().GetSize()[0];

    FloatImageType::IndexType index[8];

    index[0][2] = std::max(0, std::min(static_cast<int>(input_depth-1),  static_cast<int>(floor(z_transformed))));
    index[0][1] = std::max(0, std::min(static_cast<int>(input_height-1), static_cast<int>(floor(y_transformed))));
    index[0][0] = std::max(0, std::min(static_cast<int>(input_width-1),  static_cast<int>(floor(x_transformed))));

    index[1][2] = std::max(0, std::min(static_cast<int>(input_depth-1),  static_cast<int>(ceil( z_transformed))));
    index[1][1] = std::max(0, std::min(static_cast<int>(input_height-1), static_cast<int>(floor(y_transformed))));
    index[1][0] = std::max(0, std::min(static_cast<int>(input_width-1),  static_cast<int>(floor(x_transformed))));

    index[2][2] = std::max(0, std::min(static_cast<int>(input_depth-1),  static_cast<int>(floor(z_transformed))));
    index[2][1] = std::max(0, std::min(static_cast<int>(input_height-1), static_cast<int>(ceil( y_transformed))));
    index[2][0] = std::max(0, std::min(static_cast<int>(input_width-1),  static_cast<int>(floor(x_transformed))));

    index[3][2] = std::max(0, std::min(static_cast<int>(input_depth-1),  static_cast<int>(ceil( z_transformed))));
    index[3][1] = std::max(0, std::min(static_cast<int>(input_height-1), static_cast<int>(ceil( y_transformed))));
    index[3][0] = std::max(0, std::min(static_cast<int>(input_width-1),  static_cast<int>(floor(x_transformed))));

    index[4][2] = std::max(0, std::min(static_cast<int>(input_depth-1),  static_cast<int>(floor(z_transformed))));
    index[4][1] = std::max(0, std::min(static_cast<int>(input_height-1), static_cast<int>(floor(y_transformed))));
    index[4][0] = std::max(0, std::min(static_cast<int>(input_width-1),  static_cast<int>(ceil( x_transformed))));

    index[5][2] = std::max(0, std::min(static_cast<int>(input_depth-1),  static_cast<int>(ceil( z_transformed))));
    index[5][1] = std::max(0, std::min(static_cast<int>(input_height-1), static_cast<int>(floor(y_transformed))));
    index[5][0] = std::max(0, std::min(static_cast<int>(input_width-1),  static_cast<int>(ceil( x_transformed))));

    index[6][2] = std::max(0, std::min(static_cast<int>(input_depth-1),  static_cast<int>(floor(z_transformed))));
    index[6][1] = std::max(0, std::min(static_cast<int>(input_height-1), static_cast<int>(ceil( y_transformed))));
    index[6][0] = std::max(0, std::min(static_cast<int>(input_width-1),  static_cast<int>(ceil( x_transformed))));

    index[7][2] = std::max(0, std::min(static_cast<int>(input_depth-1),  static_cast<int>(ceil( z_transformed))));
    index[7][1] = std::max(0, std::min(static_cast<int>(input_height-1), static_cast<int>(ceil( y_transformed))));
    index[7][0] = std::max(0, std::min(static_cast<int>(input_width-1),  static_cast<int>(ceil( x_transformed))));

    const float i1 = inputImage->GetPixel(index[0]) * (1-zd) + inputImage->GetPixel(index[1]) * zd;
    const float i2 = inputImage->GetPixel(index[2]) * (1-zd) + inputImage->GetPixel(index[3]) * zd;
    const float j1 = inputImage->GetPixel(index[4]) * (1-zd) + inputImage->GetPixel(index[5]) * zd;
    const float j2 = inputImage->GetPixel(index[6]) * (1-zd) + inputImage->GetPixel(index[7]) * zd;

    const float w1 = i1 * (1-yd) + i2 * yd;
    const float w2 = j1 * (1-yd) + j2 * yd;

    // finally return the calculated value
    return (w1 * (1-xd) + w2 * xd);
  }

  /** \brief Calculates eigenvalues and corresponding eigenvectors of a symmetric 3x3 matrix
    *
    * Eigenvalues are sorted (absolute value); eigenvectors are normalized (to length = 1)
    * NOTE: Function originates from the public domain Java Matrix library JAMA.
    *
    * \param[out] V is a 3x3 matrix of the 3 eigenvectors,
    *               where the vector V[0][0], V[1][0], V[2][0] corresponds to d[0] (smallest absolute eigenvalue)
    *               norm of each vector is 1
    * \param[out] d is an array containing the 3 eigenvalues sorted by their absolute values
    *               fabs(d[0]) <= fabs(d[1]) <= fabs(d[2])
    */
  static void calculateEigenDecompositionSymmetric3x3( float A[3][3], float V[3][3], float d[3] )
  {
      float e[3];
      float da[3];
      float dt, dat;
      float vet[3];
      int i, j;
      for (i = 0; i < 3; i++) {
          for (j = 0; j < 3; j++) {
              V[i][j] = A[i][j];
          }
      }
      tred2(V, d, e);
      tql2(V, d, e);

      /* Sort the eigen values and vectors by abs eigen value */
      da[0]=fabsf(d[0]); da[1]=fabsf(d[1]); da[2]=fabsf(d[2]);
      if((da[0]>=da[1])&&(da[0]>da[2]))
      {
          dt=d[2];   dat=da[2];    vet[0]=V[0][2];    vet[1]=V[1][2];    vet[2]=V[2][2];
          d[2]=d[0]; da[2]=da[0];  V[0][2] = V[0][0]; V[1][2] = V[1][0]; V[2][2] = V[2][0];
          d[0]=dt;   da[0]=dat;    V[0][0] = vet[0];  V[1][0] = vet[1];  V[2][0] = vet[2];
      }
      else if((da[1]>=da[0])&&(da[1]>da[2]))
      {
          dt=d[2];   dat=da[2];    vet[0]=V[0][2];    vet[1]=V[1][2];    vet[2]=V[2][2];
          d[2]=d[1]; da[2]=da[1];  V[0][2] = V[0][1]; V[1][2] = V[1][1]; V[2][2] = V[2][1];
          d[1]=dt;   da[1]=dat;    V[0][1] = vet[0];  V[1][1] = vet[1];  V[2][1] = vet[2];
      }
      if(da[0]>da[1])
      {
          dt=d[1];   dat=da[1];    vet[0]=V[0][1];    vet[1]=V[1][1];    vet[2]=V[2][1];
          d[1]=d[0]; da[1]=da[0];  V[0][1] = V[0][0]; V[1][1] = V[1][0]; V[2][1] = V[2][0];
          d[0]=dt;   da[0]=dat;    V[0][0] = vet[0];  V[1][0] = vet[1];  V[2][0] = vet[2];
      }
  }

private:


  /* Eigen decomposition code for symmetric 3x3 matrices, copied from the public
   * domain Java Matrix library JAMA. */
  // little helper
  static float hypot2(float x, float y) { return sqrtf(x*x+y*y); }

  /* Symmetric Householder reduction to tridiagonal form. */
  static void tred2(float V[3][3], float d[3], float e[3])
  {
  /*  This is derived from the Algol procedures tred2 by */
  /*  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for */
  /*  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding */
  /*  Fortran subroutine in EISPACK. */
      int i, j, k;
      float scale;
      float f, g, h;
      float hh;
      for (j = 0; j < 3; j++) {d[j] = V[3-1][j]; }

      /* Householder reduction to tridiagonal form. */

      for (i = 3-1; i > 0; i--) {
          /* Scale to avoid under/overflow. */
          scale = 0.0;
          h = 0.0;
          for (k = 0; k < i; k++) { scale = scale + fabsf(d[k]); }
          if (scale == 0.0) {
              e[i] = d[i-1];
              for (j = 0; j < i; j++) { d[j] = V[i-1][j]; V[i][j] = 0.0;  V[j][i] = 0.0; }
          } else {

              /* Generate Householder vector. */

              for (k = 0; k < i; k++) { d[k] /= scale; h += d[k] * d[k]; }
              f = d[i-1];
              g = sqrtf(h);
              if (f > 0) { g = -g; }
              e[i] = scale * g;
              h = h - f * g;
              d[i-1] = f - g;
              for (j = 0; j < i; j++) { e[j] = 0.0; }

              /* Apply similarity transformation to remaining columns. */

              for (j = 0; j < i; j++) {
                  f = d[j];
                  V[j][i] = f;
                  g = e[j] + V[j][j] * f;
                  for (k = j+1; k <= i-1; k++) { g += V[k][j] * d[k]; e[k] += V[k][j] * f; }
                  e[j] = g;
              }
              f = 0.0;
              for (j = 0; j < i; j++) { e[j] /= h; f += e[j] * d[j]; }
              hh = f / (h + h);
              for (j = 0; j < i; j++) { e[j] -= hh * d[j]; }
              for (j = 0; j < i; j++) {
                  f = d[j]; g = e[j];
                  for (k = j; k <= i-1; k++) { V[k][j] -= (f * e[k] + g * d[k]); }
                  d[j] = V[i-1][j];
                  V[i][j] = 0.0;
              }
          }
          d[i] = h;
      }

      /* Accumulate transformations. */

      for (i = 0; i < 3-1; i++) {
          V[3-1][i] = V[i][i];
          V[i][i] = 1.0;
          h = d[i+1];
          if (h != 0.0) {
              for (k = 0; k <= i; k++) { d[k] = V[k][i+1] / h;}
              for (j = 0; j <= i; j++) {
                  g = 0.0;
                  for (k = 0; k <= i; k++) { g += V[k][i+1] * V[k][j]; }
                  for (k = 0; k <= i; k++) { V[k][j] -= g * d[k]; }
              }
          }
          for (k = 0; k <= i; k++) { V[k][i+1] = 0.0;}
      }
      for (j = 0; j < 3; j++) { d[j] = V[3-1][j]; V[3-1][j] = 0.0; }
      V[3-1][3-1] = 1.0;
      e[0] = 0.0;
  }

  /* Symmetric tridiagonal QL algorithm. */
  static void tql2( float V[3][3], float d[3], float e[3] )
  {
  /*  This is derived from the Algol procedures tql2, by */
  /*  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for */
  /*  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding */
  /*  Fortran subroutine in EISPACK. */

      int i, j, k, l, m;
      float f;
      float tst1;
      float eps;
      int iter;
      float g, p, r;
      float dl1, h, c, c2, c3, el1, s, s2;

      for (i = 1; i < 3; i++) { e[i-1] = e[i]; }
      e[3-1] = 0.0;

      f = 0.0;
      tst1 = 0.0;
      eps = 0.0000019209f;  //1.9209e-07 => machine epsilon => take e-06 to be on the safe side
      for (l = 0; l < 3; l++) {

          /* Find small subdiagonal element */

          tst1 = fmaxf(tst1, fabsf(d[l]) + fabsf(e[l]));
          m = l;
          while (m < 3) {
              if (fabsf(e[m]) <= eps*tst1) { break; }
              m++;
          }

          /* If m == l, d[l] is an eigenvalue, */
          /* otherwise, iterate. */

          if (m > l) {
              iter = 0;
              do {
                  iter = iter + 1;  /* (Could check iteration count here.) */
                  /* Compute implicit shift */
                  g = d[l];
                  p = (d[l+1] - g) / (2.0f * e[l]);
                  r = hypot2(p, 1.0f);
                  if (p < 0) { r = -r; }
                  d[l] = e[l] / (p + r);
                  d[l+1] = e[l] * (p + r);
                  dl1 = d[l+1];
                  h = g - d[l];
                  for (i = l+2; i < 3; i++) { d[i] -= h; }
                  f = f + h;
                  /* Implicit QL transformation. */
                  p = d[m]; c = 1.0; c2 = c; c3 = c;
                  el1 = e[l+1]; s = 0.0; s2 = 0.0;
                  for (i = m-1; i >= l; i--) {
                      c3 = c2;
                      c2 = c;
                      s2 = s;
                      g = c * e[i];
                      h = c * p;
                      r = hypot2(p, e[i]);
                      e[i+1] = s * r;
                      s = e[i] / r;
                      c = p / r;
                      p = c * d[i] - s * g;
                      d[i+1] = h + s * (c * g + s * d[i]);
                      /* Accumulate transformation. */
                      for (k = 0; k < 3; k++) {
                          h = V[k][i+1];
                          V[k][i+1] = s * V[k][i] + c * h;
                          V[k][i] = c * V[k][i] - s * h;
                      }
                  }
                  p = -s * s2 * c3 * el1 * e[l] / dl1;
                  e[l] = s * p;
                  d[l] = c * p;

                  /* Check for convergence. */
              } while (fabsf(e[l]) > eps*tst1);
          }
          d[l] = d[l] + f;
          e[l] = 0.0;
      }

      /* Sort eigenvalues and corresponding vectors. */
      for (i = 0; i < 3-1; i++) {
          k = i;
          p = d[i];
          for (j = i+1; j < 3; j++) {
              if (d[j] < p) {
                  k = j;
                  p = d[j];
              }
          }
          if (k != i) {
              d[k] = d[i];
              d[i] = p;
              for (j = 0; j < 3; j++) {
                  p = V[j][i];
                  V[j][i] = V[j][k];
                  V[j][k] = p;
              }
          }
      }
  }

};

#endif // NUMERICSHELPER_H
