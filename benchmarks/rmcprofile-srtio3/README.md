# Native RMCProfile7 SrTiO3 tutorial benchmark

This benchmark uses the official RMCProfile tutorial Exercise 4 rather than a
reconstructed publication example. It compares zero-move neutron F(Q), total
G(r), and all six partial RDFs from the identical 2,880-atom 5 K SrTiO3
configuration.

The RMCProfile7b.35 distribution still ships the exercise in legacy `.rmc6f`
format. `scripts/run_native.py` performs the format conversion with the bundled
RMCCreate utility, restricts real space to the minimum-image limit of 15.6 A,
and removes Bragg and distance-window blocks because they are outside this
forward-model comparison. It then runs the official binary at zero moves.

```bash
./fetch_reference.sh
python3 scripts/run_native.py \
  --archive reference/RMCProfile_V7b35_Linux_ARM64.tgz
cargo build --release
python3 scripts/run_comparison.py \
  --native-run-dir results/native/run
```

Docker Desktop is required because the official archive supplies a Linux ARM64
binary. The 49,500-atom `02big` configuration currently fails inside
RMCProfile7b.35 box initialization; the official 2,880-atom `01small` case runs
successfully and permits comparison through 15 A.

## Result

The pointwise neutron F(Q) comparison contains 1,963 points and gives an RMS
difference of 1.03% of the native RMCProfile curve's dynamic range, with a
maximum difference of 4.99%. Histogram peaks fall on different point labels:
RMCProfile labels each 0.02 A bin by its upper edge, while rsmith uses the bin
centre. For that reason the real-space acceptance comparison rebins both curves
to 0.1 A before pointwise scoring. The rebinned total G(r) RMS difference is
2.41% of its range; the six partial-RDF RMS differences span 0.056--4.86%.

These results demonstrate native forward parity on an official RMCProfile
tutorial configuration. The raw 0.02 A pointwise RDF errors are also retained
in `results/rsmith/summary.json` so the bin-origin difference remains visible.
Frozen regression envelopes are recorded in `expected/native-forward.toml`.
