#include <cstdint>
extern "C" void edge_correlations(const double* z,int n,int t,const int64_t* ptr,const int32_t* idx,double* corr){
 #pragma omp parallel for schedule(dynamic,32)
 for(int i=0;i<n;++i){
  const double* x=z+int64_t(i)*t;
  for(int64_t e=ptr[i];e<ptr[i+1];++e){
   const double* y=z+int64_t(idx[e])*t;double dot=0;
   #pragma omp simd reduction(+:dot)
   for(int j=0;j<t;++j)dot+=x[j]*y[j];
   corr[e]=dot/t;
  }
 }
}
extern "C" void same_network_sums(const int32_t* labels,int maps,int n,int k,const int64_t* ptr,const int32_t* idx,const double* corr,const double* area,double* sums,double* weights){
 #pragma omp parallel for schedule(dynamic,1)
 for(int m=0;m<maps;++m){
  const int32_t* lab=labels+int64_t(m)*n;double* s=sums+int64_t(m)*k;double* w=weights+int64_t(m)*k;
  for(int i=0;i<n;++i){int c=lab[i]-1;if(c<0||c>=k)continue;
   for(int64_t e=ptr[i];e<ptr[i+1];++e){int j=idx[e];if(lab[j]!=lab[i])continue;double a=area[i]*area[j];s[c]+=a*corr[e];w[c]+=a;}
  }
 }
}
