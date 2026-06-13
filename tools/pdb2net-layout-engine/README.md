# PDB2Net Headless Layout Engine

This is a small dependency-free Java CLI that calculates deterministic
Prefuse-like force-directed node positions for PDB2Net headless CX2 exports.

It is not a vendored copy of Cytoscape or Prefuse. It implements the same broad
force-directed model semantics used by Prefuse-style layouts:

- N-body repulsion
- spring force for edges
- damping
- deterministic initialization
- component-wise layout and compact packing

## Build

Install a JDK, then run:

```bash
cd tools/pdb2net-layout-engine
./build.sh
```

The jar is written to:

```bash
tools/pdb2net-layout-engine/build/pdb2net-layout-engine.jar
```

Use this path as `layout_engine_path` with:

```json
{
  "layout_mode": "prefuse_headless",
  "layout_engine_path": "tools/pdb2net-layout-engine/build/pdb2net-layout-engine.jar"
}
```

## CLI

```bash
java -jar build/pdb2net-layout-engine.jar \
  --input /path/to/layout_job.json \
  --output /path/to/positions.json
```

Input:

```json
{
  "network_title": "Example",
  "nodes": [{"id": "A"}, {"id": "B"}],
  "edges": [{"source": "A", "target": "B"}],
  "layout": {
    "algorithm": "prefuse_force_directed",
    "iterations": 100,
    "spring_coefficient": 0.0001,
    "spring_length": 50,
    "node_mass": 3,
    "use_edge_weights": false,
    "deterministic": true
  }
}
```

Output:

```json
{
  "positions": {
    "A": {"x": -80.0, "y": 0.0},
    "B": {"x": 80.0, "y": 0.0}
  }
}
```
