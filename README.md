## UnivSketch

### Introduction

Measuring flow cardinality, also known as flow spread, in large and high-speed data streams is essential for a variety of practical applications. The complex patterns within the data stream continually inspire a growing array of monitoring tasks concerning flow cardinality, including per-flow cardinality estimation, k-persistent spread estimation, super-spreader detection, and heavy cardinality changer detection, which underscores a critical need for simultaneous monitoring of multiple statistics to offer comprehensive analysis. However, it remains an unresolved challenge to deal with the diverse requirements of these tasks on the storage and retrieval of cardinality data while optimizing computational and memory efficiency. In this paper, we propose UnivSketch, the first universal measurement for measuring flow cardinality tailored for broad applicability and optimized performance across various monitoring tasks. UnivSketch features a novel virtual multi-resolution bitmap (VMRB) structure for universal encoding of flow-element data and a heavy part for tracking flows with significant cardinalities, thereby supporting the tasks of interest. Furthermore, we devise a virtual martingale mechanism that utilizes the encoded data of VMRB to deliver real-time and accurate cardinality estimates to the heavy part, serving as the basis for deciding which flows to track. We have theoretically established the error bounds, asymptotic unbiasedness, and the property of constant-time insertion for UnivSketch. Experimental results on four real-world datasets demonstrate UnivSketch's capability to concurrently monitor multiple flow cardinality metrics, surpassing baseline methods in terms of estimation accuracy and F1-Score across four distinct tasks while also achieving a notable 1.86x~13.67x improvement on insertion throughput.

### About this repo

The core **UnivSketch** structure is implemented in **/algorithms/univ**.

Other baseline methods are also implemented in **/algorithms**.

The dataset files are placed under the **/data**.

The general modules are placed under the **/utils**.

The main function is placed under the **/simulation**.

### Requirements

- g++ (gcc-version >= 13.1.0)

### Dataset

Please note that the dataset needs to have its IP addresses converted to integers beforehand. Each line should consist of two unsigned integers separated by a space: the first number representing the source IP address and the second number representing the destination IP address.

Additionally, due to the large size of the original dataset, the data we provided in **/data** only consists of unique integers pairs after deduplication.

### How to build

You can use the following commands to build and run.

```
cd simulation
g++ ./main.cpp -o main -std=c++17 -I ../../UnivSketch -mavx2 -O3
./main 2048
```

> The num '2048' after './main' specifies the allocated memory size (2048KB) .
