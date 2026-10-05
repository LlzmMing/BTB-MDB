# BTB-MDB

BTB is a C++ solver for the maximum defective biclique problem in bipartite graphs. It searches for a solution with at most `k` missing edges between its two vertex sets and at least `q` vertices on each side, maximizing the number of existing edges between the selected sets.

## Build

Run these commands from the project root (the directory containing `CMakeLists.txt`):

```bash
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -- -j4
cd ..
```

The executable is generated at `bin/run` in the project root.


## Input format

The input is a plain-text bipartite edge list:

```text
m n_left n_right
u_1 v_1
u_2 v_2
...
u_m v_m
```

- `m` is the number of edges.
- `n_left` and `n_right` are the numbers of vertices in the two partitions.
- Each pair `u v` connects left vertex `u` to right vertex `v`.
- Vertex IDs are zero-based and numbered independently in each partition: `0 <= u < n_left` and `0 <= v < n_right`.
- Provide exactly `m` edge pairs, without duplicates or comment lines. End the file with a newline.

For example, create a small dataset from the project root:

```bash
mkdir -p examples/toy
cat > examples/toy/edges.txt <<'DATA'
5 3 2
0 0
0 1
1 0
1 1
2 0
DATA
```

## Run

Pass an input file:

```bash
./bin/run -d examples/toy/edges.txt -k 1 -q 2
```

Alternatively, pass a directory containing `edges.txt`:

```bash
./bin/run -d examples/toy -k 1 -q 2
```

For this example, all three left vertices and both right vertices contain five existing edges and one missing edge, so the optimal edge count is `5`.

To request four OpenMP threads:

```bash
./bin/run -d examples/toy/edges.txt -k 1 -q 2 -j 4
```

Show command-line help:

```bash
./bin/run --help
```

### Parameters

| Option | Meaning | Default |
| --- | --- | --- |
| `-d`, `--data` | Input file, or a directory containing `edges.txt`. | Required |
| `-k`, `--key` | Maximum total number of missing edges between the selected vertex sets. Use a nonnegative integer. | Required |
| `-q`, `--lb` | Minimum number of selected vertices on each side. Use a nonnegative integer. | `0` |


## Output

Progress and timing information are printed to standard error. Each completed run appends a tab-separated row to `result.txt` in the **current working directory**, without a header:

```text
instance  input_edges  left_vertices  right_vertices  k  q  completed  optimal_edges  branches  pruned  seconds
```

`completed` is `1` for a completed run. `optimal_edges` is the reported objective value, `branches` and `pruned` are search counters, and `seconds` is elapsed solver time. For an input named `edges.txt`, the instance name is taken from its parent directory; otherwise, the input filename is used.

When the process receives `SIGTERM`, `SIGINT`, or `SIGALRM`, it appends a shorter row containing the first seven fields, with `completed` set to `0`, and exits with status `124`. The program does not impose a time limit itself. On Linux, an external limit can be applied with:

```bash
timeout 600s ./bin/run -d examples/toy/edges.txt -k 1 -q 2
```
