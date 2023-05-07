#include <cstdio>
#include <cstdlib>
#include <vector>

__global__ void bucket_sort(int *bucket, int *key, int n, int range ){
  for (int i=0; i<range; i++) {
    bucket[i] = 0;
  }
  for (int i=0; i<n; i++) {
    bucket[key[i]]++;
  }
  for (int i=0, j=0; i<range; i++) {
    for (; bucket[i]>0; bucket[i]--) {
      key[j++] = i;
    }
  }
}



int main() {
  int n = 50;
  int range = 5;
  std::vector<int> key(n);
  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]);
  }
  printf("\n");

  int* deviceKey;
  int* deviceBucket;

  cudaMalloc((void**)&deviceKey, n * sizeof(int));
  cudaMalloc((void**)&deviceBucket, range * sizeof(int));

  cudaMemcpy(deviceKey, key.data(), n * sizeof(int), cudaMemcpyHostToDevice);

  bucket_sort<<<1, 1>>>(deviceBucket, deviceKey, n, range);

  cudaMemcpy(key.data(), deviceKey, n * sizeof(int), cudaMemcpyDeviceToHost);

  cudaFree(deviceBucket);
  cudaFree(deviceKey);

  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");
}
