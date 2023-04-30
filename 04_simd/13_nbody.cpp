#include <cstdio>
#include <cstdlib>
#include <cmath>

int main() {
  const int N = 8;
  float x[N], y[N], m[N], fx[N], fy[N];
  for(int i=0; i<N; i++) {
    x[i] = drand48();
    y[i] = drand48();
    m[i] = drand48();
    fx[i] = fy[i] = 0;
  }
  __m256 xvec = _mm256_load_ps(x);
  __m256 yvec = _mm256_load_ps(y);
  __m256 mvec = _mm256_load_ps(m);
  __m256 fxvec = _mm256_load_ps(fx);
  __m256 fyvec = _mm256_load_ps(fy);

  __m256 zero = _mm256_setzero_ps();
  __m256 one = _mm256_set1_ps(1);

  for(int i=0; i<N; i++) {
    __m256 xi_vec = _mm256_set1_ps(x[i]);
    __m256 yi_vec = _mm256_set1_ps(y[i]);

    __m256 rx = _mm256_sub_ps(xvec,xi_vec);
    __m256 ry = _mm256_sub_ps(yvec,yi_vec);

    __m256 r2 = _mm256_add_ps(_mm256_mul_ps(rx, rx), _mm256_mul_ps(ry, ry));
    __m256 r = _mm256_sqrt_ps(r2);

    __m256 mask = _mm256_cmp_ps(r, zero, _CMP_GT_OQ);

    r2 =  _mm256_blendv_ps(one, r2 , mask);
    r = _mm256_blendv_ps(one,r,mask);

    __m256 dfx = _mm256_div_ps(_mm256_mul_ps(rx,mvec),_mm256_mul_ps(r2,r));
    __m256 dfy = _mm256_div_ps(_mm256_mul_ps(ry,mvec),_mm256_mul_ps(r2,r));


    _mm256_store_ps(fx,dfx);
    _mm256_store_ps(fy,dfy);

    for(int j=0; j<N; j++){
        fx[i] -= fx[j];
        fy[i] -= fx[j];
    }
    printf("%d %g %g\n",i,fx[i],fy[i]);
  }
}
