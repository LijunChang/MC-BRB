# MC-BRB
Computing Maximum Clique in Large Sparse Graphs

# Compile
```sh
$ make clean
$ make
```

# Usage
```sh
./MC-BRB {alg} {graph_directory}
```
alg is chosen from "MC-DD, MC-EGO, MC-BRB, verify"

For data format, please see https://github.com/LijunChang/Cohesive_subgraph_book/tree/master/datasets

# Note
If you have any question, please contact me at ljchang.au@gmail.com.

If you are using this code, please kindly cite the following paper.

Lijun Chang: "Efficient Maximum Clique Computation over Large Sparse Graphs". KDD 2019: 529-538

```
@inproceedings{DBLP:conf/kdd/Chang19,
  author    = {Lijun Chang},
  title     = {Efficient Maximum Clique Computation over Large Sparse Graphs},
  booktitle = {Proceedings of the 25th {ACM} {SIGKDD} International Conference on
               Knowledge Discovery {\&} Data Mining, {KDD} 2019, Anchorage, AK,
               USA, August 4-8, 2019.},
  pages     = {529--538},
  year      = {2019},
  crossref  = {DBLP:conf/kdd/2019},
  url       = {https://doi.org/10.1145/3292500.3330986},
  doi       = {10.1145/3292500.3330986},
  timestamp = {Sun, 28 Jul 2019 15:40:22 +0200},
  biburl    = {https://dblp.org/rec/bib/conf/kdd/Chang19},
  bibsource = {dblp computer science bibliography, https://dblp.org}
}
```
