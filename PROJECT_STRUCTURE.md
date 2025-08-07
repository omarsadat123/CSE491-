# Dense-PCE Project Structure

## Overview
This repository contains implementations of various pseudo-clique enumeration algorithms, including Dense-PCE, PCE, FPCE, and EBBkC. The project focuses on finding dense subgraphs (pseudo-cliques) in large-scale networks with different optimization strategies.

## Repository Structure

```
CSE491-/
├── README.md                                    # Main project README
├── WSL-Ubuntu-cmds.txt                         # WSL Ubuntu commands reference
└── Dense_PCE/                                  # Main project directory
    └── Dense-PCE-main/                         # Core implementation directory
        ├── dense-pce.cpp                       # Main Dense-PCE implementation (C++)
        ├── dense-pce                           # Compiled Dense-PCE executable
        ├── dense-pce.exe                       # Windows executable
        ├── readme.md                           # Dense-PCE usage instructions
        ├── run.sh                              # Shell script for running experiments
        ├── run_both.sh                         # Script to run both PCE and Dense-PCE
        ├── transgrh.pl                         # Graph format conversion utility
        ├── Graph_check.py                      # Graph preprocessing utility
        ├── generate_synth_graphs.ipynb         # Synthetic graph generation notebook
        ├── parse output log.ipynb              # Output analysis notebook
        │
        ├── Algorithms/                         # Algorithm implementations
        │   ├── pce12/                          # Original PCE implementation
        │   │   ├── pce.c                       # Main PCE algorithm
        │   │   ├── pce                         # Compiled PCE executable
        │   │   ├── problem.c/h                 # Problem definition
        │   │   ├── sgraph.c/h                  # Graph data structures
        │   │   ├── queue.c/h                   # Queue implementation
        │   │   ├── itemset.c/h                 # Itemset handling
        │   │   ├── vec.c/h                     # Vector operations
        │   │   ├── aheap.c/h                   # Heap implementation
        │   │   ├── alist.c/h                   # Adjacency list
        │   │   ├── base.c/h                    # Base utilities
        │   │   ├── undo.c/h                    # Undo functionality
        │   │   ├── stdlib2.c/h                 # Extended standard library
        │   │   ├── makefile                    # Build configuration
        │   │   ├── readme.txt                  # PCE documentation
        │   │   └── *.pl                        # Perl utilities for graph processing
        │   │
        │   ├── FPCE/                           # Fast Pseudo-Clique Enumeration
        │   │   ├── code/                       # FPCE source code
        │   │   │   ├── FPCE/                   # Main FPCE implementation
        │   │   │   ├── PCE/                    # PCE variant
        │   │   │   ├── ODES/                   # Dense subgraph algorithm
        │   │   │   ├── SIGMOD24-MQCE-main/     # FastQC implementation
        │   │   │   └── cluster_one-1.0.jar     # ClusterONE tool
        │   │   ├── real-synth-graph-analysis/  # Graph analysis experiments
        │   │   ├── e-coli-analysis/            # E. coli data analysis
        │   │   ├── README.md                   # FPCE documentation
        │   │   └── *.pl                        # Graph conversion utilities
        │   │
        │   └── EBBkC/                          # Efficient k-clique listing
        │       ├── src/                        # EBBkC source code
        │       │   ├── truss/                  # Truss decomposition
        │       │   ├── dependencies/           # External libraries
        │       │   │   ├── cub/                # CUDA utilities
        │       │   │   ├── libpopcnt/          # Population count library
        │       │   │   ├── moderngpu/          # Modern GPU utilities
        │       │   │   ├── sparsehash-yche/    # Sparse hash implementation
        │       │   │   └── sparsepp/           # Sparse++ library
        │       │   └── util/                   # Utility functions
        │       ├── dataset/                    # Test datasets
        │       ├── Cohesive_subgraph_book/     # Graph algorithms reference
        │       ├── README.md                   # EBBkC documentation
        │       └── EBBkC_TR.pdf                # Technical report
        │
        ├── Data/                               # Graph datasets
        │   ├── testGraphs/                     # Test graph files
        │   │   └── fpce_graph/                 # FPCE test graphs
        │   ├── synth-graphs-1000/              # Synthetic graphs (1000 nodes)
        │   ├── graph.grh                       # Large graph file
        │   ├── bn-fly-drosophila_medulla_1.grh # Drosophila brain network
        │   └── cs3                             # Graph file
        │
        ├── Results/                            # Experimental results
        │   ├── csv_output_multi_nodes_multi_desnity_100_1000_001_30.csv
        │   ├── csv_output_multi_nodes_multi_desnity_100_500_001_35.csv
        │   ├── output_both_log_1000.txt
        │   ├── output_cs3_log_1000.txt
        │   ├── output_log.txt
        │   ├── output_log_1000.txt
        │   ├── output_log_10000.txt
        │   └── nohup.out
        │
        └── Documentation/                      # Documentation files
            ├── 08pce_journal.pdf              # PCE journal paper
            └── fpce.pdf                        # FPCE paper
```

## Algorithm Implementations

### 1. Dense-PCE (Main Implementation)
- **File**: `dense-pce.cpp`
- **Language**: C++
- **Purpose**: Efficient pseudo-clique enumeration with density-based optimization
- **Key Features**:
  - Uses degree-based tracking system
  - Implements sophisticated branching strategy
  - Supports configurable theta threshold
  - Includes Turan filtering capabilities
  - Core decomposition for seed clique generation

### 2. PCE (Original Implementation)
- **Directory**: `pce12/`
- **Language**: C
- **Purpose**: Original pseudo-clique enumeration algorithm
- **Key Features**:
  - Itemset-based approach
  - Queue-based processing
  - Graph format conversion utilities
  - Comprehensive graph data structures

### 3. FPCE (Fast Pseudo-Clique Enumeration)
- **Directory**: `FPCE/`
- **Language**: C
- **Purpose**: Optimized pseudo-clique enumeration with edge bounding
- **Key Features**:
  - Edge-bound optimization
  - Multiple graph format support
  - Real and synthetic graph analysis
  - E. coli data analysis capabilities

### 4. EBBkC (Efficient k-clique Listing)
- **Directory**: `EBBkC/`
- **Language**: C++
- **Purpose**: k-clique enumeration with edge-oriented branching
- **Key Features**:
  - CUDA support for GPU acceleration
  - Truss decomposition
  - Early termination strategies
  - Parallel processing capabilities

## Build Instructions

### Dense-PCE
```bash
cd Dense_PCE/Dense-PCE-main
g++ -std=c++11 -O3 dense-pce.cpp -o dense-pce
```

### PCE
```bash
cd Dense_PCE/Dense-PCE-main/pce12
make
mv pce ..
```

### FPCE
```bash
cd Dense_PCE/Dense-PCE-main/FPCE/code/FPCE
make
```

### EBBkC
```bash
cd Dense_PCE/Dense-PCE-main/EBBkC/src
mkdir build && cd build
cmake ..
make
```

## Usage Examples

### Dense-PCE
```bash
./dense-pce synth-graphs-1000/scale_free_graph_m_10.grh --maximum 300 --theta 1.0
```

### PCE
```bash
./pce C -l 1 -u 100 bn-fly-drosophila_medulla_1.grh 0.9
```

### FPCE
```bash
./code/FPCE/fpce M -l 10 real-graphs/bio-grid-human.grh 0.9
```

### EBBkC
```bash
./BBkC e dataset/facebook 20 3
```

## Graph Formats Supported

1. **EDGELIST**: Space-separated node pairs
2. **GRH**: Adjacency list format (0 to n-1 node IDs)
3. **Binary**: Optimized binary format (b_adj.bin, b_degree.bin)

## Key Features

### Graph Processing
- Multiple format support and conversion utilities
- Preprocessing for self-loops and duplicate edges
- Node ID normalization

### Algorithm Optimizations
- Degree-based tracking systems
- Early termination strategies
- Edge bounding techniques
- Parallel processing support

### Analysis Capabilities
- Real-world graph analysis
- Synthetic graph generation
- Performance benchmarking
- Memory usage tracking

### Experimental Framework
- Automated experiment scripts
- Result logging and analysis
- CSV output generation
- Jupyter notebook integration

## Dependencies

### System Requirements
- **Operating System**: Linux (Ubuntu 20.04 recommended)
- **Compiler**: GCC 9.4.0 or higher
- **Java**: OpenJDK 11.0.22 (for ClusterONE)
- **CMake**: Version 3.16+ (for EBBkC)
- **CUDA**: For GPU acceleration (optional)

### External Libraries
- **CUB**: CUDA utilities
- **ModernGPU**: GPU programming utilities
- **SparseHash**: Memory-efficient hash tables
- **Sparse++**: Sparse data structures
- **libpopcnt**: Population count operations

## Research Context

This repository implements and compares multiple approaches to pseudo-clique enumeration:

1. **Dense-PCE**: Novel density-based approach with sophisticated branching
2. **PCE**: Original pseudo-clique enumeration algorithm
3. **FPCE**: Fast variant with edge bounding optimizations
4. **EBBkC**: k-clique enumeration with early termination

The implementations support research in:
- Graph mining and analysis
- Dense subgraph discovery
- Network community detection
- Biological network analysis
- Social network analysis

## Current Status

- **Dense-PCE**: Fully implemented and functional
- **PCE**: Original implementation, stable
- **FPCE**: Complete with experimental framework
- **EBBkC**: GPU-accelerated k-clique enumeration
- **Documentation**: Comprehensive usage instructions
- **Experiments**: Automated testing and benchmarking scripts 