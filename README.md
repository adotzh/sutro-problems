# sutro-problems

A collection of small, self-contained problems used as benchmarks and
illustrations for the Sutro Group's energy-efficient learning research.

## Structure

Each problem lives in **its own folder**, with **its own README** describing:

- **The problem** — what is being learned, the input/output format, and the dataset
- **Reference algorithm** — the original method as published / proposed
- **Reproduction code** — runnable scripts (Python or otherwise)
- **Reproduction results** — expected accuracy, training time, and any plots
- **Notes** — variants, follow-ups, open questions

To add a new problem, create a new directory and add a `README.md` at its root.

## Problems

- [`wip-boltzmann-shifter/`](wip-boltzmann-shifter/) — the shift-direction
  inference task from Hinton & Sejnowski's Boltzmann-machine chapter
  (PDP Vol 1, 1986). _(work in progress)_
- [`energy-efficient-matmul/`](energy-efficient-matmul/) — cheapest matrix
  multiplication under a simplified Dally explicit-communication cost
  model (data movement priced; arithmetic free). Baselines: 4×4 naive
  1,316; 16×16 naive 340,704; 16×16 tiled 133,783.
