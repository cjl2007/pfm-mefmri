// Exact local edge corrections. No dense connectivity matrix is constructed.
#include <cstdint>
extern "C" void ridge_local(const double* z, int nc, int nt,
 const int64_t* ptr, const int32_t* idx, const double* a, int k,
 double* correction, double* grams) {
 #pragma omp parallel for schedule(dynamic,32)
 for (int v=0; v<nc; ++v) {
  const double* x=z+int64_t(v)*nt;
  double* c=correction+int64_t(v)*k;
  double* g=grams+int64_t(v)*k*k;
  for (int64_t e=ptr[v]; e<ptr[v+1]; ++e) {
   int j=idx[e]; const double* y=z+int64_t(j)*nt; double dot=0;
   #pragma omp simd reduction(+:dot)
   for (int t=0; t<nt; ++t) dot+=x[t]*y[t];
   dot/=nt; const double* row=a+int64_t(j)*k;
   for (int u=0; u<k; ++u) {
    c[u]+=dot*row[u];
    for(int w=0; w<=u; ++w) g[u*k+w]+=row[u]*row[w];
   }
  }
  for(int u=0;u<k;++u) for(int w=0;w<u;++w) g[w*k+u]=g[u*k+w];
 }
}
