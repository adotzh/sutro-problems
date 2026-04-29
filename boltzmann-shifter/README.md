# boltzmann-shifter

The **shift-direction inference** task from Hinton & Sejnowski, *"Learning and
Relearning in Boltzmann Machines"* (Chapter 7 of Rumelhart, McClelland & PDP
Research Group, *Parallel Distributed Processing*, Vol 1, MIT Press, 1986).

## Problem

Given two rings of `N` binary units `V1` and `V2`, where `V2` is a copy of
`V1` shifted (with wraparound) by one of `{-1, 0, +1}` positions, infer
which shift was applied.

- **Input**: `V1` (N bits), `V2` (N bits)
- **Output**: one of three classes — `left`, `none`, `right`
- **Training set**: full enumeration of `2^N` patterns × `3` shifts
- **For N=8**: 768 training cases

The interesting property of the problem is that no single bit of `V2` is
sufficient to identify the shift — the network must discover **multiplicative
position-pair detectors** between `V1` and `V2`.

## Files

| File | Purpose |
|---|---|
| `shifter.py` | Fully-connected Boltzmann machine with simulated-annealing Gibbs sampling. Faithful to the 1986 procedure but very slow to converge in finite wall-clock time. |
| `shifter_rbm.py` | Restricted Boltzmann Machine trained with CD-1 (Hinton 2002). Same fundamental learning rule (positive-phase minus negative-phase statistics), but with the efficient RBM sampling structure. **This is the working version.** |
| `visualize_shifter.py` | Generates per-unit receptive-field panels, full weight heatmap, hidden-activation maps, and confusion matrix. |
| `figure3.py` | Renders the trained network's hidden units in the Hinton-diagram style (white/black squares sized by `|w|`) — same layout convention as Figure 3 of the original chapter. |
| `viz/` | Output PNGs from the runs below. |

## Running

```bash
python3 shifter_rbm.py --N 8 --hidden 80 --epochs 200 --lr 0.03 --momentum 0.7 --batch 32
```

Training takes ~110 seconds on a laptop. Final accuracy on the full 768-case
test set: **86.98%**.

To regenerate all visualization outputs:

```bash
python3 visualize_shifter.py --N 8 --hidden 64 --epochs 200 --outdir viz
python3 figure3.py --N 8 --epochs 300 --outdir viz
```

## Results

| Metric | Value |
|---|---|
| Final accuracy (full 768 cases) | 86.98% |
| Per-class: `left (-1)` | 92.6% (237/256) |
| Per-class: `none (0)` | 84.0% (215/256) |
| Per-class: `right (+1)` | 86.3% (221/256) |
| Training time | ~110 sec |
| Hyperparameters | hidden=80, lr=0.03, momentum=0.7, batch=32, 200 CD-1 epochs |

The hidden-unit population partitions cleanly by preferred shift class
(roughly 1/3 each for left/none/right, mirroring the 3-class structure of
the task). Many units learn the expected position-pair detectors: a unit
preferring "shift left" tends to have correlated weights at `V1[i]` and
`V2[i-1 mod N]`.

## Deviations from the 1986 procedure

1. **Architecture** — RBM (visible-hidden only) instead of fully-connected.
   Same learning rule applied to a sparser graph; permits exact one-step
   Gibbs sampling.
2. **Sampling** — CD-1 (Hinton 2002) instead of simulated annealing.
3. **Hidden units** — 80 here (24 in the original, with most of those
   "doing very little" per the original analysis). The fully-connected
   `shifter.py` does run with 24 hidden units to match the original Figure 3
   layout.
4. **Hardware** — modern laptop, ~minutes; the original ran on a VAX with
   substantially longer training time.

## Open questions / next experiments

- Does the receptive-field structure that emerges from CD-1 quantitatively
  match the structure reported for the 1986 simulated-annealing run?
- How does the metric we use to evaluate algorithms (FLOPs, ByteDMD,
  data-movement cost) rank these two training procedures? The annealing
  schedule is dominated by sampling work; CD-1 is dominated by matrix
  multiplies.
- Can the shifter task be solved with a direct contrastive objective (no
  generative pre-training), and how does its sample efficiency compare?
