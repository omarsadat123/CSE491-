# Dense-PCE

### compile PCE
```bash
cd pce12
make
mv pce ..
cd ..
```

### compile Dense-PCE

```bash
g++ -O3 dense-pce.cpp -o dense-pce 
```

### generate graph

generate graph using generate_synth_graphs.ipynb notebook.


### Run PCE
```bash
./pce  C -l <lower limit of pseudo-clique size> -u <upper limit of pseudo-clique size> <path-to-grh-file> <theta>

# example
./pce  C -l 1 -u 100 bn-fly-drosophila_medulla_1.grh .9
```

### Run Dense-PCE
```bash
./dense-pce <path-to-grh-file>  --maximum <maximum-possible-size-of-pseudo-clique> --theta <theta>

# example
./dense-pce synth-graphs-1000/scale_free_graph_m_10.grh --maximum 300 --theta 1.0
```