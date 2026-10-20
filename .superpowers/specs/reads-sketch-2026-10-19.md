# Spec: `repgenr reads-sketch` and `reads-screen`, stream, sketch and judge sequencing runs (2026-10-19)

Status: implemented. `reads-sketch` in PR #35, the variant layout in #36,
`reads-screen` with `assemble --accession-file` in #37, the Nextflow
`READS_SKETCH` and `READS_SCREEN` processes in #38 (all merged 2026-10-19).
This file records the design as shipped; the user documentation is in
`docs/usage.md` (sections "Sketching runs without downloading them",
"Screening the sketches before any download", "Reads in Nextflow") and
`docs/output.md` (`reads-sketch`, `reads-screen`).

## Purpose

Sketch the reads of ENA/SRA runs without downloading them. The reads are
streamed from the archive with `fastq-dump` and piped into one `sourmash
sketch dna` call, so no FASTQ file touches the disk and a run can be cut
after its first N spots. The result is one reads sketch per run in the form
the pipeline already uses (k=21, 31, 51, scaled=1000, with abundance), ready
for a `sourmash gather` species screen, duplicate detection and surveys over
many runs. `reads-screen` applies the screen of `assemble --screen-reads` to
those sketches, so runs can be left out before any download.

Decisions taken with the user on 2026-10-19:

- `reads` finds the runs; `reads-sketch` works on a workdir's `reads.tsv`,
  all runs or a selection of them, and records the result in the workdir.
- Separate subcommands, not part of the assemble download path. `assemble`
  keeps fetching whole files and sketching them beside the assembler.
- The spot cap is optional. Without it the whole run is streamed.
- Each cap is a variant directory of its own, so sketches at several depths
  coexist and the number of spots is visible in the layout.
- Streaming by `fastq-dump` (sra-tools), not by HTTP range requests on the
  ENA mirror. `fasterq-dump` has no spot cap (checked, sra-tools 3.4.1).

## Commands

The flow is three commands: `reads` selects the runs of a taxon or a study
and writes `reads.tsv`; `reads-sketch` sketches those runs, all of them or a
subset, without the download; `reads-screen` judges the sketches and writes
the accepted runs to a file `assemble --accession-file` takes.

```
repgenr reads -wd WD -ts "Francisella tularensis" --platform illumina
repgenr reads-sketch -wd WD --max-spots 200000
repgenr reads-sketch -wd WD --accession SRR1 --accession SRR2        # a subset
repgenr reads-screen -wd WD --max-spots 200000 --gtdb-sketch DB --gtdb-lineages LIN
repgenr assemble -wd WD --accession-file WD/reads_screen/spots_200000/passed.txt
```

### `reads-sketch`

| Option | Meaning |
|---|---|
| `-wd, --workdir` | A workdir whose `reads` stage has completed. Required unless `-o` is given. |
| `--accession ACC` | Restrict to these runs of `reads.tsv`, repeatable. A run not in `reads.tsv` is refused by name. |
| `--accession-file FILE` | The same, one accession per line, `#` comments allowed. Combines with `--accession`. |
| `--platform` | `any`, `illumina`, `ont`, `pacbio`: restrict by the platform column of `reads.tsv`. Needs `-wd`. |
| `--max-spots N` | Stream only spots 1..N of each run (`fastq-dump -X N`). A paired-end spot gives two reads. Default: the whole run. Names the variant. |
| `--jobs J` | Runs streamed at once. Default `min(4, usable_cpus())`. Not part of the resume fingerprint. |
| `--force` | Stream again the selected runs of the current variant. |
| `-o, --out DIR` | Standalone form: sketch into `DIR/<variant>/` without a workdir; nothing is recorded. Takes bare `--accession` runs or `--reads-tsv`. |
| `--reads-tsv FILE` | Standalone form: the runs of a `reads.tsv` (the Nextflow task's input). |
| `--versions-out FILE` | Standalone form: the tool versions as a YAML fragment (Nextflow). |

Without a selection every run of `reads.tsv` is sketched. A call adds to
its variant: runs sketched before with the same cap are kept, so a subset
can be sketched first and the rest later.

### `reads-screen`

| Option | Meaning |
|---|---|
| `-wd, --workdir` | A workdir with `reads.tsv` and `reads_sketches/`. Required unless the standalone form is used. |
| `--accession`, `--accession-file` | As for `reads-sketch`. |
| `--max-spots N` | The variant to screen (`reads_sketches/spots_N/`). Without it `all/` when present, else the only variant there; with several, required. |
| `--min-fraction F` | Minimum abundance-weighted fraction of the reads in the top GTDB species (default 0.3). |
| `--dup-ani A` | ANI estimate from which a run contained in an earlier run of the same biosample is a duplicate (default 0.999). |
| `--gtdb-sketch`, `--gtdb-lineages` | The GTDB sketch and its lineages CSV, or `REPGENR_GTDB_SKETCH` and `REPGENR_GTDB_LINEAGES`. |
| `--jobs J` | Gathers at once (default `min(4, usable_cpus())`, about 0.6 GB each). |
| `--reads-tsv`, `--sketches-dir DIR`, `-o OUT` | Standalone form: screen the runs of a `reads.tsv` from the `-o DIR` of `reads-sketch` into `OUT/<variant>/`; nothing is recorded. |
| `--versions-out FILE` | Standalone form: the tool versions (Nextflow). |

### `assemble --accession ACC`, `--accession-file FILE`

The runs of `reads.tsv` to assemble; the other pending runs are deferred
as under `--max-runs`, finished runs are kept, an accession outside
`reads.tsv` is refused, and a file that lists no accessions is refused
(an empty `passed.txt` means the screen excused every run).

Panel: the entry panel, `reads-sketch` and `reads-screen` after `reads`.

## Behaviour

One pipeline per run:

```
fastq-dump -Z --split-spot --skip-technical [-X N] ACC \
  | sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund --name ACC -o <tmp> -
```

- `-Z` writes to stdout, `--split-spot` emits both mates of a pair as
  separate records, and `--skip-technical` drops barcode and adapter reads.
- The parameter string is `READS_SKETCH_PARAMS` from `core/sketches.py`, so
  the file is interchangeable with the `reads.sig.zip` the assemble stage
  writes. Abundance is kept because the screen's gather weights by it.
- The signature is written to a hidden temporary file beside the target and
  renamed into place, so an interrupted run leaves no partial sketch.
- `fastq-dump` reports `Read N spots` and `Written N spots` on stderr and
  sourmash `calculated 1 signature for N sequences`; both go into the index.
- Each run streams in a disposable `.home.<run>.*` working directory inside
  the variant directory, so the `.ncbi` settings `fastq-dump` writes under
  `HOME` (the task cwd under Docker) do not land in the output.
- Failure attribution: `fastq-dump` is blamed only when its stderr tail
  holds `err: ` or `quit with error code`; otherwise sourmash when it exited
  non-zero (it exiting first makes `fastq-dump` fail on the closed pipe
  with an empty stderr).
- Retries (2026-10-20 audit): a failure attributed to `fastq-dump` whose
  message is not permanent (not found, withdrawn, protected, no data), or a
  pipeline past its time limit, is tried again, three attempts with pauses
  of about 30 s and 120 s (a quarter of jitter). With a cap the per-run
  limit is 300 s plus `max_spots / 200` s, or `REPGENR_SUBPROCESS_TIMEOUT`
  when lower; without one only the global limit applies. Worst case per
  run at 200000 spots: three attempts of 1300 s plus the pauses, about
  68 min, before its slot is freed for good.
- A stream with a spot count, no error line, and no sequence for sourmash
  (`no sequences found`, exit 0, no file: no spots, or technical reads
  only) is the status `empty`: settled at that
  cap (`is_settled`), no deliverable, not streamed again without
  `--force`; `is_current` still means a readable sketch, so `reads-screen`
  excuses the run as unsketched.
- One call per variant directory: `<variant>/.lock` is held with `flock`
  for the call and the index is read under it; a second call is refused
  with a message; a file system without locks gives a warning.
- Both images are pulled once before the pool (`containers.ensure_image`);
  a Docker exit 125 with the engine's own error line stops the call as one
  error (`containers.engine_failure_line`).
- The sourmash half runs with `RAYON_NUM_THREADS=1` through
  `run_pipeline(tool_envs=...)`; `tool_env` is the route by which a
  variable reaches the tool itself in every backend (Docker `-e`).
- The per-run working directory is `.home.<run>.<random>`, its cleanup
  best effort.

Measured 2026-10-19 on this Mac: 5000 spots of a MiSeq run streamed and
sketched in 12 to 14 s. Under `--container docker` each tool runs in its own
image with the host pipe between them. Audit 2026-10-20 (20 runs at 200000
spots, Docker, 4 jobs): about 10 s of fixed cost per run at NCBI, about
10 s of stream, sourmash 39 s emulated (the BioContainer is amd64-only) or
14 s on the host with the default thread pool and 7 s with one thread; 8
streams against 4 gained 15 percent; one probe stalled for 160 s.

### Output

```
WD/
  reads.tsv                               from the reads stage, unchanged
  reads_sketches/
    all/                                  the whole run
      <run>.sig.zip
      reads_sketches.tsv                  the index of this variant
    spots_<N>/                            --max-spots N
      <run>.sig.zip
      reads_sketches.tsv
  reads_screen/
    <variant>/
      reads_screen.tsv                    one row per run of reads.tsv
      passed.txt                          the accepted runs, reads.tsv order
      <run>/screen.json                   one decision with its settings
      <run>/screen/                       the gather and tax genome output
```

With `-o DIR` the same layout goes to `DIR/<variant>/`.

`reads_sketches.tsv`, one row per run of `reads.tsv` in its order, rewritten
after every finished run: `run_accession`, `status` (`ok`, `failed`,
`pending`), `max_spots`, `spots_read`, `sequences`, `sketch_file`,
`sketch_params`, `seconds`, `reason`.

`reads_screen.tsv` has the columns of `screen_reads.tsv`: `run_accession`,
`top_match`, `top_genus`, `fraction`, `duplicate_of`, `decision` (`pass`,
`excused`), `reason` (`taxon_mismatch`, `low_match_fraction`,
`duplicate_isolate`, or `screen_failed` for a run without a sketch or one
the tools failed on).

### Stage records and resume

- `reads_sketch`: inputs `reads.tsv`, the variant's index and the accession
  file; deliverables the index and the sketch of every selected run, so a
  failed run is a missing deliverable and a repeat streams it again; in
  `_REDIGEST_AFTER_RUN` because the index is its own input. Params record
  the variant of the call and the runs each variant holds; `status` lists
  the variants. A call where every selected run failed raises before the
  record is written, so the stage stays incomplete.
- `reads_screen`: inputs `reads.tsv`, the variant's index and the accession
  file; deliverables the table, `passed.txt` and the decision of every
  selected run with a sketch (a `screen_failed: gather` run keeps no
  decision, so it is retried). A decision made with the same settings on
  the same sketch (by mtime) is reused; a changed threshold or a newer
  sketch screens again from the kept sketch. The duplicate check is redone
  on every call in `reads.tsv` order over every run with a decision made
  with the same settings, selected or not, so several calls give the
  result of one call. A run decided with other settings and not selected
  keeps its row but leaves `passed.txt`, with a warning. Rows of runs no
  longer in `reads.tsv` are dropped.

### Tools and containers

- `core/sratools.py`: `ToolCapabilities` `sra-tools` with
  `fastq-dump >= 3.0`, conda `bioconda::sra-tools`, the image
  `ncbi/sra-tools:3.4.1` (NCBI's multi-architecture image; the amd64-only
  BioContainer segfaulted under emulation on Apple Silicon); `sra-tools>=3.0` in
  `envs/core.yml`.
- `core/process.run_pipeline` and `core/containers.run_pipeline` join two
  tools by a pipe, each half with its own image, environment and `--bin-dir`,
  both registered as live processes; when the consumer fails the producer's
  process group is stopped and both containers are removed.
- `reads-screen` reuses `assemble_screen.ScreenGate` unchanged; the settings
  builder `assemble_screen.resolve_settings` is shared with `assemble`.

### What this is not

- Not a download. No FASTQ is kept; `assemble` fetches whole files.
- Not a random sample. The first N spots of an Illumina run come from the
  first tiles, those of an ONT run from the start of the run; a few hundred
  thousand spots are a reasonable starting point, not a tested threshold.
- With a cap, the duplicate check compares the starts of two runs rather
  than the runs; the estimate is less certain than from whole runs.
- `assemble --screen-reads` does not read these decisions: its sketch is of
  the whole run.

## Nextflow

`READS_SKETCH`, `READS_SKETCH_MERGE` and `READS_SCREEN` in
`modules/local/dataflow/`, wired into `ACQUIRE_READS` between `READS_SELECT`
and the per-run `READS_ASSEMBLE` scatter. The sketch runs one task per
chunk of `params.reads_sketch_chunk` runs (default 25; the task takes
`reads.tsv` and the chunk's accessions through `--accession-file`), and the
merge joins the chunk directories into one `sketches/` with one index per
variant in `reads.tsv` order (shell and awk, no tool). The screen is one
task per batch. `params.reads_sketch` runs the sketch
task (`params.reads_sketch_spots` is the cap) and publishes
`reads/reads_sketches/<variant>/`; `params.reads_screen` implies it, runs
the screen with `params.gtdb_sketch`, `params.gtdb_lineages` and
`params.reads_screen_args`, publishes `reads/reads_screen/<variant>/`, and
only the runs of `passed.txt` become assembly tasks. A screen that passes no
run stops the pipeline; `reads_screen` without the GTDB parameters is
refused before any task (the pipeline does not read the `REPGENR_GTDB_*`
environment variables). Both tasks run the standalone command forms. Stub
nf-tests cover both processes and the reads data flow with and without the
screen.

## Resolved open points

1. Command names: `reads-sketch`, `reads-screen`.
2. Without `fastq-dump` the command is refused (preflight); the ENA mirror
   answers HTTP range requests and a gzip prefix decodes to whole reads
   (verified 2026-10-19), but that second code path was not built.
3. `reads-sketch` and `reads-screen` select by accession list and platform
   only; `reads` applies the size filters when it writes `reads.tsv`.
