# 🧬 OptJV3

**OptJV3** is a thread-parallel genetic algorithm optimizer for discovering high-performing **JARVIS3** compression parameter sets. It treats JARVIS3 as a black-box compressor and searches over global mixer parameters, context models, and repeat models to minimize compressed output size or related compression objectives.

OptJV3 is intended for reproducible compression experiments on DNA or other sequence datasets where the parameter space is too large for exhaustive search.

---

## ✨ Highlights

* 🧬 Genetic algorithm search over the JARVIS3 command-line parameter space.
* ⚡ Thread-parallel evaluation of candidate configurations.
* 🧠 Configurable context-model search through `-cm` genes.
* 🔁 Configurable repeat-model search through `-rm` genes.
* 🌐 Search or fix global JARVIS3 parameters: `-hs`, `-lr`, and `-sd`.
* 🎯 Optional restart from a known-good parameter string.
* ✅ Optional decompression verification.
* 📈 Per-generation CSV history output.
* 🏆 Best-solution report with reusable JARVIS3 command line.
* 🧪 Support for release, debug, and PGO builds through the included Makefile.

---

## 📁 Repository layout

```text
.
├── Makefile                       # Build system for OptJV3
├── optjv3.c                       # OptJV3 source code
├── OptJV3                         # Compiled optimizer binary
├── JARVIS3                        # JARVIS3 executable
├── cy_sample                      # Example input sample
├── RunOptJV3.sh                   # Standard optimization run
├── RunOptJV3FromSpecificPoint.sh  # Restart run from a known parameter set
├── best.txt                       # Best solution written by OptJV3
├── history.csv                    # Per-generation optimization history
└── LICENSE                        # GPL v3 license text
```

The examples below assume that `OptJV3`, `JARVIS3`, and `cy_sample` are located in the repository root.

---

## 🛠️ Requirements

OptJV3 is written in **C11** and targets POSIX-compatible systems.

Required:

* `gcc`
* `make`
* POSIX shell
* POSIX threads
* math library, linked with `-lm`
* compiled `JARVIS3` executable

The default Makefile uses:

```makefile
CC       := gcc
TARGET   := OptJV3
SRC      := optjv3.c
STD      := -std=c11
WARN     := -Wall -Wextra
THREADS  := -pthread
OPT      := -O3 -march=native -flto=auto
LIBS     := -lm
```

---

## 🚀 Quick start

Build the release binary:

```bash
make
```

or explicitly:

```bash
make release
```

Make sure both executables are available and executable:

```bash
chmod +x ./OptJV3 ./JARVIS3
```

Run the standard optimization example:

```bash
./RunOptJV3.sh
```

---

## 🧱 Building

The project is built using the included **Makefile**.

### 🏁 Release build

```bash
make release
```

This builds:

```text
OptJV3
```

using:

```text
-O3 -march=native -flto=auto -std=c11 -Wall -Wextra -pthread -lm
```

### 🐞 Debug build

```bash
make debug
```

The debug target cleans previous objects and rebuilds with:

```text
-O0 -g3 -std=c11 -Wall -Wextra -pthread
```

Use this build for debugging with tools such as `gdb` or `valgrind`.

### 🔁 Rebuild

```bash
make rebuild
```

Equivalent to:

```bash
make clean
make release
```

### 🧹 Clean

Remove the binary and object files:

```bash
make clean
```

Remove build artifacts and PGO profile data:

```bash
make distclean
```

### 🧪 Profile-guided optimization

The Makefile includes a two-stage **PGO** workflow.

First, build an instrumented binary:

```bash
make pgo-gen
```

Then run representative workloads, for example:

```bash
./RunOptJV3.sh
```

After profile data have been collected in `profdir`, rebuild using the profiles:

```bash
make pgo-use
```

Remove PGO data:

```bash
make pgo-clean
```

---

## 📜 Example scripts

### ▶️ `RunOptJV3.sh`

This script performs a compact optimization run with a population of 32 individuals and 10 generations.

```bash
#!/usr/bin/env bash
set -euo pipefail

./OptJV3 \
  --jarvis ./JARVIS3 \
  --input cy_sample \
  --population 32 \
  --generations 10 \
  --threads 4 \
  --elite 4 \
  --tournament 3 \
  --crossover 0.85 \
  --mutation 0.12 \
  --toggle 0.08 \
  --objective bytes \
  --max-cmodels 2 \
  --min-cmodels 2 \
  --max-rmodels 1 \
  --min-rmodels 1 \
  --global-bounds "hs=8:32,lr=0.0:0.15,seed=1:1000" \
  --cm-bounds "ctx=1:13,den=1:200,ir=0:2,gamma=0.01:0.99,edits=0:3,eden=1:50,eir=0:1,egamma=0.01:0.99" \
  --rm-bounds "nr=1:20,ctx=11:13,beta=0.01:0.99,limit=1:20,gamma=0.01:0.99,ir=0:1,weight=0.01:0.99,cache=1:4" \
  --best-out best.txt \
  --history-out history.csv
```

### 🎯 `RunOptJV3FromSpecificPoint.sh`

This script seeds the population from a known parameter configuration and searches nearby configurations. It is useful after a previous run has found a strong solution.

```bash
#!/usr/bin/env bash
set -euo pipefail

./OptJV3 \
  --jarvis ./JARVIS3 \
  --input cy_sample \
  --population 100 \
  --generations 100 \
  --threads 8 \
  --elite 4 \
  --tournament 3 \
  --crossover 0.80 \
  --mutation 0.12 \
  --toggle 0.08 \
  --objective bytes \
  --max-cmodels 3 \
  --min-cmodels 2 \
  --max-rmodels 2 \
  --min-rmodels 1 \
  --global-bounds "hs=8:32,lr=0.0:0.0,seed=1:1000" \
  --cm-bounds "ctx=1:13,den=1:500,ir=0:2,gamma=0.01:0.99,edits=0:3,eden=1:50,eir=0:1,egamma=0.01:0.99" \
  --rm-bounds "nr=1:20,ctx=11:13,beta=0.01:0.99,limit=1:20,gamma=0.01:0.99,ir=0:1,weight=0.01:0.99,cache=1:4" \
  --best-out best.txt \
  --history-out history.csv \
  --restart-from "./JARVIS3 -hs 21 -lr 0.000 -sd 17 -cm 4:437:1:0.864/0:32:1:0.347 -cm 13:97:1:0.872/1:12:0:0.880 -rm 17:11:0.933:3:0.631:0:0.059:2 -rm 18:11:0.738:3:0.428:1:0.010:1"
```

`--restart-from` accepts a JARVIS3-style parameter string. OptJV3 parses `-hs`, `-lr`, `-sd`, `-cm`, and `-rm`. A leading executable token such as `./JARVIS3` is harmless and ignored.

---

## ⚙️ How OptJV3 works

OptJV3 maintains a population of candidate JARVIS3 configurations. Each candidate is converted into a JARVIS3 command line, executed on the input file, and scored according to the selected objective.

At each generation:

1. Candidate parameter strings are generated or inherited.
2. JARVIS3 compresses the input file with each candidate.
3. Fitness is computed from compressed size, bits per symbol, or size plus runtime.
4. Candidates are sorted by fitness.
5. The best candidates survive through elitism.
6. New candidates are created using tournament selection, crossover, mutation, and model toggling.
7. The generation best is printed and optionally appended to `history.csv`.

Lower fitness is better.

---

## 🧩 Candidate representation

Each candidate contains global parameters, context-model genes, and repeat-model genes.

### 🌐 Global parameters

| Parameter | JARVIS3 flag | Meaning                                        |
| --------- | -----------: | ---------------------------------------------- |
| `hs`      |        `-hs` | Hidden-size / neural mixer capacity parameter. |
| `lr`      |        `-lr` | Learning rate.                                 |
| `seed`    |        `-sd` | JARVIS3 seed.                                  |

### 🧠 Context model gene

Context models use the following format:

```text
ctx:den:ir:gamma/edits:eden:eir:egamma
```

Example:

```text
-cm 13:97:1:0.872/1:12:0:0.880
```

| Field    | Meaning                                    |
| -------- | ------------------------------------------ |
| `ctx`    | Context depth.                             |
| `den`    | Denominator / smoothing-related parameter. |
| `ir`     | Inverted-repeat mode.                      |
| `gamma`  | Decay or mixing parameter.                 |
| `edits`  | Number of edit contexts.                   |
| `eden`   | Edit denominator parameter.                |
| `eir`    | Edit inverted-repeat mode.                 |
| `egamma` | Edit-model decay or mixing parameter.      |

### 🔁 Repeat model gene

Repeat models use the following format:

```text
nr:ctx:beta:limit:gamma:ir:weight:cache
```

Example:

```text
-rm 17:11:0.933:3:0.631:0:0.059:2
```

| Field    | Meaning                                     |
| -------- | ------------------------------------------- |
| `nr`     | Number/range parameter for repeat modeling. |
| `ctx`    | Repeat-model context depth.                 |
| `beta`   | Repeat-model beta parameter.                |
| `limit`  | Repeat limit.                               |
| `gamma`  | Repeat-model decay or mixing parameter.     |
| `ir`     | Inverted-repeat mode.                       |
| `weight` | Repeat-model weight.                        |
| `cache`  | Cache level.                                |

The repeat-model order is fixed throughout OptJV3:

```text
nr:ctx:beta:limit:gamma:ir:weight:cache
```

---

## 🧾 Command-line options

### ✅ Required options

| Option          | Description                       |
| --------------- | --------------------------------- |
| `--jarvis PATH` | Path to the JARVIS3 executable.   |
| `--input FILE`  | Input file used for optimization. |

### 🧬 Genetic algorithm options

| Option               |    Default | Description                                                        |
| -------------------- | ---------: | ------------------------------------------------------------------ |
| `--population N`     |       `48` | Number of individuals in the population.                           |
| `--generations N`    |       `30` | Number of generations to run.                                      |
| `--threads N`        |        `4` | Number of parallel worker threads.                                 |
| `--elite N`          |        `4` | Number of best candidates copied unchanged to the next generation. |
| `--tournament N`     |        `3` | Tournament size used for parent selection.                         |
| `--crossover X`      |     `0.85` | Probability of crossover when creating offspring.                  |
| `--mutation X`       |     `0.12` | Per-field mutation probability.                                    |
| `--toggle X`         |     `0.08` | Probability of enabling/disabling a model slot.                    |
| `--blend X`          |     `0.35` | Blend-alpha value for real-valued crossover.                       |
| `--seed U64`         | time-based | Master RNG seed for OptJV3.                                        |
| `--restart-from STR` |   disabled | Seed the initial population from a known parameter string.         |

### 🎚️ Objective options

| Option             | Values                       | Description                                      |
| ------------------ | ---------------------------- | ------------------------------------------------ |
| `--objective NAME` | `bytes`, `bps`, `bytes+time` | Fitness objective. Lower is better.              |
| `--time-weight X`  | floating point               | Time penalty weight used only with `bytes+time`. |

Objective definitions:

* `bytes`: minimize compressed file size in bytes.
* `bps`: minimize bits per input symbol.
* `bytes+time`: minimize `compressed_bytes + time_weight * elapsed_seconds`.

### 🔢 Model-count options

| Option            | Default | Description                                                       |
| ----------------- | ------: | ----------------------------------------------------------------- |
| `--max-cmodels N` |     `4` | Maximum number of context-model slots available to the optimizer. |
| `--min-cmodels N` |     `0` | Minimum number of enabled context models.                         |
| `--max-rmodels N` |     `2` | Maximum number of repeat-model slots available to the optimizer.  |
| `--min-rmodels N` |     `0` | Minimum number of enabled repeat models.                          |

Set `min=max` to force an exact number of active models.

Example:

```bash
--max-cmodels 2 --min-cmodels 2 --max-rmodels 1 --min-rmodels 1
```

This forces exactly two context models and one repeat model.

### 🌐 Global parameter control

| Option             | Default | Description                                 |                            |
| ------------------ | ------: | ------------------------------------------- | -------------------------- |
| `--optimize-hs 0   |      1` | `1`                                         | Whether to optimize `-hs`. |
| `--optimize-lr 0   |      1` | `1`                                         | Whether to optimize `-lr`. |
| `--optimize-seed 0 |      1` | `0`                                         | Whether to optimize `-sd`. |
| `--fixed-hs N`     |    `42` | Fixed `-hs` value when `--optimize-hs 0`.   |                            |
| `--fixed-lr X`     |  `0.03` | Fixed `-lr` value when `--optimize-lr 0`.   |                            |
| `--fixed-seed N`   |    `17` | Fixed `-sd` value when `--optimize-seed 0`. |                            |

### 📏 Bounds

Bounds are passed as comma-separated `key=min:max` lists.

#### 🌐 Global bounds

```bash
--global-bounds "hs=8:32,lr=0.0:0.15,seed=1:1000"
```

Supported keys:

```text
hs, lr, seed
```

#### 🧠 Context-model bounds

```bash
--cm-bounds "ctx=1:13,den=1:200,ir=0:2,gamma=0.01:0.99,edits=0:3,eden=1:50,eir=0:1,egamma=0.01:0.99"
```

Supported keys:

```text
ctx, den, ir, gamma, edits, eden, eir, egamma
```

#### 🔁 Repeat-model bounds

```bash
--rm-bounds "nr=1:20,ctx=11:13,beta=0.01:0.99,limit=1:20,gamma=0.01:0.99,ir=0:1,weight=0.01:0.99,cache=1:4"
```

Supported keys:

```text
nr, ctx, beta, limit, gamma, ir, weight, cache
```

---

## 📤 Output files

### 🏆 `best.txt`

Written when `--best-out best.txt` is used.

Contains the best solution found, including:

* fitness
* compressed bytes
* bits per symbol
* elapsed time
* validity flag
* objective name
* JARVIS3 parameter string
* full executable command

Example structure:

```text
fitness=1234.00 bytes=1234 bps=1.2345 elapsed=0.42 valid=1 objective=bytes
params: -hs 21 -lr 0.000 -sd 17 -cm ... -rm ...
command: ./JARVIS3 -hs 21 -lr 0.000 -sd 17 -cm ... -rm ... cy_sample
```

### 📈 `history.csv`

Written when `--history-out history.csv` is used.

Contains one row per generation:

```csv
generation,fitness,bytes,bps,elapsed_s,valid,objective,params
```

This file can be used to plot convergence curves, compare runs, or inspect previous experiments.

---

## 🏃 Runtime behavior

For each candidate, OptJV3 creates an isolated temporary directory under `--workdir`, runs JARVIS3 compression, reads the compressed output size, and removes temporary files unless `--keep-temps` is enabled.

Default temporary directory:

```text
/tmp
```

Use a custom workspace:

```bash
--workdir ./tmp-optjv3
```

Keep candidate temporary directories and logs:

```bash
--keep-temps
```

Suppress warning messages for failed candidates:

```bash
--quiet
```

Enable decompression verification:

```bash
--verify
```

When `--verify` is enabled, OptJV3 decompresses each candidate output and checks that the decompressed file is identical to the original input. Invalid candidates receive a very large fitness penalty.

---

## 🔬 Reproducibility

OptJV3 supports reproducible genetic search through `--seed`.

Example:

```bash
./OptJV3 \
  --jarvis ./JARVIS3 \
  --input cy_sample \
  --seed 12345 \
  --population 32 \
  --generations 10
```

For strict reproducibility, keep the following fixed:

* OptJV3 source code
* Makefile build target and compiler flags
* JARVIS3 binary
* input file
* OptJV3 `--seed`
* population size
* generation count
* number of threads
* all bounds
* objective settings

Repeated candidates are cached by parameter string and are not recomputed during a run.

---

## 🧭 Recommended workflows

### 1. 🔎 Broad exploratory search

Use wider bounds and a moderate population:

```bash
./OptJV3 \
  --jarvis ./JARVIS3 \
  --input cy_sample \
  --population 64 \
  --generations 30 \
  --threads 8 \
  --objective bytes \
  --best-out best.txt \
  --history-out history.csv
```

### 2. 🎯 Focused refinement

After a good result is found, copy the `params:` line from `best.txt` and pass it to `--restart-from`.

Use tighter bounds and more generations:

```bash
./OptJV3 \
  --jarvis ./JARVIS3 \
  --input cy_sample \
  --population 100 \
  --generations 100 \
  --threads 8 \
  --objective bytes \
  --restart-from "-hs 21 -lr 0.000 -sd 17 -cm 4:437:1:0.864/0:32:1:0.347 -cm 13:97:1:0.872/1:12:0:0.880 -rm 17:11:0.933:3:0.631:0:0.059:2"
```

### 3. ⏱️ Speed-aware optimization

Use `bytes+time` when a slightly larger compressed file is acceptable if it is produced faster:

```bash
./OptJV3 \
  --jarvis ./JARVIS3 \
  --input cy_sample \
  --objective bytes+time \
  --time-weight 10.0
```

### 4. 🧪 PGO-tuned optimizer build

Use the PGO Makefile targets when repeatedly running OptJV3 on representative workloads:

```bash
make pgo-gen
./RunOptJV3.sh
make pgo-use
```

Then run the optimized binary normally:

```bash
./RunOptJV3.sh
```

---

## 📊 Interpreting results

During execution, OptJV3 prints the best candidate from each generation:

```text
=> Generation 5 of 10 best: fitness=1234.00 bytes=1234 bps=1.2345 time=0.42 valid=1
./JARVIS3 -hs 21 -lr 0.000 -sd 17 -cm ... -rm ... cy_sample
```

Important fields:

| Field     | Meaning                                               |
| --------- | ----------------------------------------------------- |
| `fitness` | Value minimized by the selected objective.            |
| `bytes`   | Compressed output size.                               |
| `bps`     | Bits per input symbol.                                |
| `time`    | Compression runtime for that candidate.               |
| `valid`   | `1` if candidate evaluation succeeded, `0` otherwise. |

At the end of the run, the final best solution is printed under:

```text
Best solution found
===================
```

---

## 🚧 Hard limits

OptJV3 includes internal hard limits to match JARVIS3 parser constraints and avoid invalid command generation.

### 🧱 Internal model-slot limits

| Constant           | Value |
| ------------------ | ----: |
| `HARD_MAX_CMODELS` |  `16` |
| `HARD_MAX_RMODELS` |   `8` |

### 🧠 Context-model parser limits

| Field             |      Range |
| ----------------- | ---------: |
| `ctx`             |    `1..14` |
| `den`             |  `1..5000` |
| `ir`              |     `0..2` |
| `edits`           |    `0..20` |
| `eden`            | `1..50000` |
| `eir`             |     `0..1` |
| `gamma`, `egamma` |   `(0, 1)` |

### 🔁 Repeat-model parser limits

| Field                     |       Range |
| ------------------------- | ----------: |
| `nr`                      | `1..100000` |
| `ctx`                     |     `1..14` |
| `limit`                   |     `1..20` |
| `ir`                      |      `0..2` |
| `cache`                   |     `1..15` |
| `beta`, `gamma`, `weight` |    `(0, 1)` |

User-provided bounds are sanitized against these limits.

---

## 🧯 Troubleshooting

### `JARVIS executable not found or not executable`

Check that the path passed to `--jarvis` exists and has execute permission:

```bash
chmod +x ./JARVIS3
```

### `Input file not readable`

Check that the input file exists and is readable:

```bash
ls -lh cy_sample
```

### Build fails with unsupported `-march=native` or LTO options

Edit the Makefile and simplify:

```makefile
OPT := -O3
```

Then rebuild:

```bash
make clean
make release
```

### Many candidates fail

Try one or more of the following:

* Narrow the search bounds.
* Use `--verify` to detect invalid decompressions explicitly.
* Use `--keep-temps` to inspect candidate logs.
* Temporarily remove `--quiet` to see warnings.
* Confirm that all generated `-cm` and `-rm` configurations are accepted by your JARVIS3 build.

### Optimization converges too early

Increase diversity:

```bash
--population 100 --mutation 0.15 --toggle 0.10
```

Or widen the parameter bounds.

### Optimization is too slow

Reduce evaluation cost:

```bash
--population 32 --generations 10 --threads 8
```

Also consider using a smaller representative input sample during exploratory search.

---

## 📌 Notes on compression benchmarking

For meaningful comparisons:

* Use the same input file across runs.
* Keep the same JARVIS3 binary.
* Record the full OptJV3 command line.
* Record the Makefile target used to build OptJV3.
* Compare final `bytes` and `bps`, not only fitness, especially when using `bytes+time`.
* Validate final candidates on independent datasets to reduce overfitting to `cy_sample`.

---

## 📚 Citation

If you use OptJV3 in academic or technical work, cite the repository and include the exact commit hash, command line, input dataset, build target, and best parameter string used for the reported result.

Example:

```text
OptJV3 genetic optimizer for JARVIS3 parameter search.
Repository: <repository-url>
Commit: <commit-hash>
Build: make release
Command: <full OptJV3 command>
Best parameters: <params line from best.txt>
```

---

## ⚖️ License

This project is licensed under the **GNU General Public License v3.0**.

You may copy, distribute, and modify this software under the terms of the GPL v3. See the [`LICENSE`](./LICENSE) file for the full license text.

Recommended SPDX identifier:

```text
GPL-3.0-only
```

---

## 🙏 Acknowledgements

OptJV3 was developed to automate parameter exploration for JARVIS3 compression experiments, especially where exhaustive search is infeasible because of a large mixed discrete/continuous parameter space.
Please, cite
```
Ferrolho, Rita, Armando J. Pinho, and Diogo Pratas. "Optimizing Genomic Data Compression with Genetic Algorithms." bioRxiv (2025): 2025-10.
```

