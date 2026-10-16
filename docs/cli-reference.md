# Command reference

Generated from the command tree by `scripts/render_cli_matrix.py`; the
matrix test keeps it in sync. Commands are grouped as in `repgenr --help`.
Global options go before the command name
(`repgenr --container docker dereplicate ...`).

A few option names differ between commands. `--tool` selects the
dereplicator on `run` and the derep commands but the SNP typer on `snptype`
(`--snptyper` on `phylo`, `run` and `phylo-build`). `run --metadata-source`
and `run --viral-source` correspond to `metadata --source` and `vmetadata
--source`. The global `--platform` is the container platform, while
`run --platform` filters the sequencing platform. On the step commands
`--outgroup-accession` takes a file that names the accession.

## Global options

| option | default | description |
|---|---|---|
| `--version` | off | Show version and exit. |
| `--container` | `none` | Run external tools in containers: none, docker, or singularity. |
| `--container-engine` |  | Engine binary override (e.g. apptainer, podman). |
| `--container-cache` |  | Directory for Singularity .sif images and their cache (large; can be external). |
| `--platform` |  | Container platform, e.g. linux/amd64 for emulated BioContainers on arm64. |
| `--wave`, `--no-wave` | off | Resolve images for multi-tool adapters via the Seqera Wave CLI. |
| `--bin-dir` |  | Put DIR first on PATH for one tool only, e.g. gubbins=ENV/bin for a satellite conda environment. Repeatable; also REPGENR_BIN_DIRS='tool=dir,tool=dir'. |
| `--force`, `-f`, `--no-force` | off | Re-run a stage even if it already completed with the same parameters. |
| `--verbose`, `-v` | off | Verbose (DEBUG) logging. |
| `--quiet`, `-q` | off | Only warnings and errors. |

## Pipeline

### run

Run the whole pipeline end to end (bacterial by default, --viral for viruses).

Some option names differ from the single commands. Here --tool selects the
dereplicator, while 'snptype --tool' selects the SNP typer (--snptyper on
phylo, run and phylo-build). --metadata-source is 'metadata --source' and --viral-source
is 'vmetadata --source'. --platform filters the sequencing platform with
--reads, while the global --platform sets the container platform.

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory (created). |
| `--viral` | off | Run the viral chain (vmetadata -> vgenome) instead of bacterial. |
| `--genomes-dir` |  | Start from local genome FASTAs (the ingest chain) instead of downloading. |
| `--from-workdir` |  | Start the ingest chain from the genome set of an earlier working directory (selection.tsv and genomes/); repeatable, and may be combined with --genomes-dir. |
| `--selection` |  | With --genomes-dir: selection.tsv naming the genomes to take (accession, taxonomy, filename, outgroup flag, quality). |
| `--outgroup` |  | With --genomes-dir or --from-workdir: the outgroup genome, a name found in those sources or a path to a FASTA file. |
| `--copy` | off | With --genomes-dir or --from-workdir: copy the files into genomes/ instead of linking. |
| `--reads` | off | Run the reads chain (reads -> assemble) from ENA/SRA sequencing runs selected by -tf/-tg/-ts or --accession-file, instead of downloading assemblies. |
| `--accession-file` |  | With --reads: file of run/sample/study accessions. |
| `--platform` | `any` | With --reads: any, illumina, ont or pacbio. |
| `--max-runs` |  | With --reads: keep at most N runs, the largest by bases. |
| `--assembler` | `auto` | Assembler: auto, flye, shovill, skesa. |
| `-d`, `--dataset` | `rep` | GTDB dataset: all or rep. |
| `-l`, `--level` |  | family, genus or species. |
| `-tf`, `--target-family` |  | Restrict the selection to this family. |
| `-tg`, `--target-genus` |  | Restrict the selection to this genus. |
| `-ts`, `--target-species` |  | Restrict the selection to this species. |
| `-r`, `--release` |  | GTDB release (tsv source). |
| `--gtdb-version` |  | GTDB table: bac120 or ar53 (tsv source). |
| `--metadata-source` | `tsv` | tsv or api. |
| `--metadata-path` |  | Use this GTDB metadata table instead of downloading. With --source tsv, -r/--release and --gtdb-version are still required. |
| `--nodownload` | off | Reuse a GTDB table already present in the workdir. |
| `--outgroup-accession` |  | Accession to fetch and set aside as the outgroup. |
| `--limit` |  | Keep at most N GTDB genomes, round-robin over species by CheckM quality. |
| `--target` |  | Virus taxon (viral). |
| `--viral-source` | `ncbi_virus` | ncbi_virus or bvbrc. |
| `--complete-only` | off | ncbi_virus: only COMPLETE sequences (viral). |
| `--host` |  | ncbi_virus: restrict to a host species (viral). |
| `--released-after` |  | ncbi_virus: MM/DD/YYYY (viral). |
| `--group-segments` | off | Group viral segments. |
| `--keep-files` | off | Keep the genome-download intermediates (bacterial chain). |
| `--tool` | `skder` | auto, drep, galah, skder, sourmash. |
| `--primary-ani` | `0.9` | Primary (pre-clustering) ANI threshold in (0, 1]. |
| `--secondary-ani` | `0.99` | Secondary (final cluster) ANI threshold in (0, 1]. |
| `--aligned-fraction` | `0.5` | Minimum aligned fraction in (0, 1] for a pair to be compared. |
| `--keeper` | `quality` | Representative choice per cluster: quality (manifest CheckM values and genome N50), gtdb (a GTDB species representative first, then quality) or tool (adapter's own). |
| `-s`, `--process-size` |  | Chunk size; when set and exceeded, two-stage chunking runs for any tool. |
| `-p`, `--num-processes` | `0` | Parallel workers, threads split across them: dereplicate's stage-1 chunk workers and, with --msa-source coregene, the genomes annotated at once. 0 = auto (~threads/4, capped by cores). |
| `--pre-primary-ani` |  | Stage-1 (intra-chunk) primary ANI; defaults to --primary-ani. |
| `--pre-secondary-ani` |  | Stage-1 (intra-chunk) secondary ANI; defaults to --secondary-ani. |
| `--reduce` | `none` | Taxonomy-aware reduction after ANI: none, species, or genus (one representative per taxon). |
| `--target-reps` | `0` | Target representative count: search --secondary-ani to land near it (0 = off; re-runs dereplication per search step). |
| `--refine` |  | Seed of a lineage to keep at finer resolution (filename, stem or accession; repeatable): every cluster holding a seed is dereplicated again at --refine-ani. |
| `--refine-ani` |  | Secondary ANI for the clusters of the --refine seeds; above --secondary-ani. |
| `--tool-arg` |  | Tool tuning as key=value (repeatable), e.g. mode=greedy. |
| `--with-snptype` | off | Run the standalone snptype stage (with --snptyper, --mask, --reference) after dereplication, so the SNP tables are produced whatever builds the tree. |
| `--treebuilder` | `iqtree` | auto, fasttree, iqtree, mashtree, raxmlng, ska2, skani, sourmash, veryfasttree. |
| `--msa-source` | `aligner` | aligner, snptype or coregene (a concatenated core-gene alignment from the annotations the annotate stage wrote). |
| `--aligner` | `progressivemauve` | cactus, progressivemauve, sibeliaz. |
| `--snptyper` | `simple` | SNP typer: parsnp, simple, ska2, snippy. |
| `--no-outgroup` | off | Do not root with an outgroup. |
| `--root` |  | How phylo roots tree/tree.nwk: outgroup (on the branch to the outgroup), midpoint, or none (as the tree builder wrote it). Default: outgroup when an outgroup is staged and --no-outgroup is not given, else none. |
| `--all-genomes` | off | Use all genomes, not only the representatives. |
| `-B`, `--bootstrap` | `0` | Bootstrap replicates. 0 turns bootstrapping off; IQ-TREE needs at least 1000 when it is on. |
| `--reference` |  | Reference genome: a filename, 'medoid' (the default: the genome nearest to all others by sketch ANI, when the genomes are sketched and sourmash is at hand) or 'first'. |
| `--aligner-arg` |  | Aligner, SNP typer, masker or tree builder tuning as key=value (repeatable), e.g. kmer=15 (sibeliaz), seed_weight=11 (progressivemauve) or model=GTR+G (iqtree, raxmlng). |
| `--mask` | `none` | Recombination masking of the SNP alignment: none, gubbins. Needs --msa-source snptype or --with-snptype. |
| `--node-basename` |  | Name internal nodes <basename><n>. Without it, internal nodes receive names derived from a hash of their descendant leaves. |
| `--naming` | `leafset` | How internal nodes are named without --node-basename: 'leafset' (a hash of the whole leaf set) or 'lca' (a hash of the smallest leaf label under each child, so most names survive genomes added with later-sorting accessions). |
| `--root-name` | `root` | Label of the top node. FlexTaxD expects 'root'. |
| `--remove-outgroup` | off | Leave the outgroup out of the taxonomy after rooting. |
| `--include-dereplicated`, `--no-include-dereplicated` | on | List redundant genomes under their representative in the taxonomy. |
| `--core-fraction` |  | Keep only the variable columns where at least this fraction of the genomes have a called base (a soft core, 0 to 1). Unset keeps every variable column and the log reports what 0.95 would keep. |
| `--coregene` | `panaroo` | Core-gene alignment builder for --msa-source coregene: panaroo, ppanggolin, scarap or orthofinder. Tuning goes through --aligner-arg. panaroo: threshold= is the identity at which it groups genes (0.98; lower it, 0.95 then 0.90, when the core comes out small on a divergent set), and also filtered=false, merge_paralogs, remove_invalid_genes, family_threshold, len_dif_percent. ppanggolin: identity= and coverage= of its protein clustering (0.8 each), partition=softcore\|core\|persistent (softcore, at --core-threshold), single_copy=false, dup_margin, mode, nb_of_partitions, translation_table. scarap: method (FH), seeds (100), max_core_genes (0, no cap), core_prefilter (the smaller of 0.90 and --core-threshold). orthofinder: search=diamond\|diamond_ultra_sens\|blastp\|mmseqs (diamond), inflation (MCL, 1.2). |
| `--annotator` | `bakta` | Annotator: bakta (default; needs its database, --bakta-db or BAKTA_DB) or prokka. |
| `--bakta-db` |  | bakta database directory (bakta_db download --output DIR --type light\|full); BAKTA_DB when unset. |
| `--core-threshold` |  | With --msa-source coregene: a gene is core when present in at least this fraction of the genomes (0.95). |
| `--clean-mode` |  | With --msa-source coregene: Panaroo's cleaning of annotation errors, strict (default), moderate or sensitive. |
| `--collapse-support` |  | Merge nodes whose support is below this fraction into their parent. |
| `--collapse-length` |  | Merge nodes whose branch is shorter than this length into their parent. |
| `--collapse-diameter` |  | Merge every clade whose leaf-to-leaf distance is within this value into one node, keeping the fewest such clades (in the tree's own units). |
| `--ranks` | off | Add a rank column to tree2tax.tsv: 'strain' on every genome, 'no rank' on internal nodes unless --rank-red names them. |
| `--rank-red` |  | With --ranks: NAME=RED bands, e.g. species=0.95,genus=0.8; an internal node takes the rank with the largest threshold not above its relative evolutionary divergence, once per lineage. Needs branch lengths and a rooted tree. |
| `-t`, `--threads` | `16` | Threads for the external tool (default 16, or the CPU limit when lower). |
| `--allow-incomplete` | off | Proceed with a warning when the input genome set is incomplete. |
| `--dry-run` | off | Print the stages and key parameters, then exit. |

### status

Show which pipeline stages have completed in a working directory.

A completed stage is listed as stale when one of its inputs changed or
one of its outputs is missing since it finished (it re-runs on its next
invocation), and as interrupted when it did not finish. A -wd that does
not exist exits 3; an existing directory without repgenr.yaml prints
which entry stage to run first and exits 0. With --json, a malformed
record still exits 3 and leaves stdout empty.

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory. |
| `--json` | off | Print one versioned JSON object on stdout instead of the text report (schema in docs/output.md); exit codes are unchanged. |

## Entry points: select and fetch genomes

### metadata

Select a taxon's genomes from GTDB (full table or the GTDB API).

The --source option is called --metadata-source in 'run'.

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory (created). |
| `-d`, `--dataset` | required | GTDB dataset: all or rep. |
| `-l`, `--level` | required | family, genus or species. |
| `--source` | `tsv` | tsv (download full table) or api (GTDB API, target only). |
| `-r`, `--release` |  | GTDB release (tsv source). |
| `--gtdb-version` |  | GTDB table: bac120 or ar53 (tsv source). |
| `-tf`, `--target-family` |  | Restrict the selection to this family. |
| `-tg`, `--target-genus` |  | Restrict the selection to this genus. |
| `-ts`, `--target-species` |  | Restrict the selection to this species. |
| `--outgroup-accession` |  | Accession to fetch and set aside as the outgroup. |
| `--metadata-path` |  | Use this GTDB metadata table instead of downloading. With --source tsv, -r/--release and --gtdb-version are still required. |
| `--nodownload` | off | Reuse a GTDB table already present in the workdir. |
| `--limit` |  | Keep at most N GTDB genomes, round-robin over species by CheckM quality. |
| `--drop-foreign` | off | Discard genomes appended from sequencing runs (assemble --append) instead of refusing to overwrite the selection that holds them. |

### genome

Download and organize genomes selected by the metadata stage.

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory. |
| `--accession-list-only` | off | Write the accession list and stop (no download). |
| `--keep-files` | off | Keep download and scratch intermediates. |
| `--sketch`, `--no-sketch` | auto | Write a sourmash sketch of each genome to sketches/ (k=21,31,51, scaled=1000). Default: when sourmash can run; --sketch requires it, --no-sketch skips it. |

### vmetadata

Retrieve viral metadata from NCBI Virus (default) or BV-BRC.

The --source option is called --viral-source in 'run'.

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory (created). |
| `--target` |  | Virus taxon/group/family. |
| `--source` | `ncbi_virus` | ncbi_virus (NCBI Virus via datasets) or bvbrc. |
| `--filter` |  | bvbrc: keep records whose header carries this tag (default: complete genome). |
| `--host` |  | ncbi_virus: restrict to a host species. |
| `--complete-only` | off | ncbi_virus: only COMPLETE sequences. |
| `--released-after` |  | ncbi_virus: MM/DD/YYYY. |
| `--list` | off | List BV-BRC targets and exit. |

### vgenome

Select and organize viral genomes (virus equivalent of genome).

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory. |
| `-tg`, `--target-genus` |  | Restrict the selection to this genus. |
| `-ts`, `--target-species` |  | Restrict the selection to this species. |
| `-tse`, `--target-serotype` |  | Restrict the selection to this serotype. |
| `-tc`, `--target-custom` |  | key:value. |
| `--length-all` | off | Disable the length window (keep every length). |
| `--length-deviation` | `10` | Half-width of the length window, in percent of its center. |
| `--length-method` | `median_of_medians` | Center of the length window: median_of_medians (default) gives one vote per species, so an over-sequenced outbreak species cannot shift the window; mean averages every record and loses that defense. |
| `--length-range` |  | e.g. 25000-35000. |
| `--discard` |  | Comma-separated header tags. |
| `--no-outgroup` | off | Do not root with an outgroup. |
| `--outgroup-accession` |  | Use this downloaded record as the outgroup instead of searching for one (an accession on the NCBI Virus path, a record id on BV-BRC). |
| `--group-segments` | off | ncbi_virus: combine an isolate's segments into one genome (segmented viruses). |
| `--outgroup-candidates-taxid-min-genomes` | `5` | Genomes a sister taxid needs to qualify as an outgroup candidate. |
| `--outgroup-treebuilder` | `mashtree` | Tree builder used for the outgroup distance matrix. Accepted: mashtree, skani. |
| `--glance` | off | Print selection and stop. |
| `--print-fasta-headers` | off | Print the headers of the selected records. |
| `--ignore-duplicates` | off | bvbrc: tolerate duplicate record ids (last wins). |
| `--keep-files` | off | Keep download and scratch intermediates. |
| `--sketch`, `--no-sketch` | auto | Write a sourmash sketch of each genome to sketches/ (k=21,31,51, scaled=1000). Default: when sourmash can run; --sketch requires it, --no-sketch skips it. |

### ingest

Populate a working directory from local genomes (no download).

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory (created). |
| `--genomes-dir` |  | Directory of genome FASTA files (.fasta, .fa, .fna, .fas, .fasta.gz, .fna.gz or .fa.gz; subdirectories are not searched) to stage under genomes/. At least one of --genomes-dir and --from-workdir is required. |
| `--from-workdir` |  | Earlier working directory whose genome set (selection.tsv and genomes/) is added, with taxonomy, quality, the GTDB representative flag and the manifest source of each genome; its outgroup is not carried over. Repeatable; may be combined with --genomes-dir. |
| `--strict-sources` | off | Refuse (exit 2, nothing written) a --from-workdir whose genome set was written by a stage that status lists as interrupted or stale, or that has no stage record; without it each such source is named in a warning and taken. |
| `--selection` |  | selection.tsv (accession, taxonomy, filename, outgroup flag, quality) naming the genomes to take from --genomes-dir; default: every FASTA under --genomes-dir, taxonomy parsed from canonical Family_genus_species_ACCESSION names. |
| `--outgroup` |  | Genome to set aside as the outgroup: a filename, stem or accession under --genomes-dir or in a --from-workdir (genomes/ or outgroup/), or a path to a FASTA file elsewhere. |
| `--copy` | off | Copy the files into genomes/ instead of symlinking them. |
| `--drop-foreign` | off | Discard genomes appended from sequencing runs (assemble --append) instead of refusing to overwrite the selection that holds them. |
| `--sketch`, `--no-sketch` | auto | Write a sourmash sketch of each genome to sketches/ (k=21,31,51, scaled=1000). Default: when sourmash can run; --sketch requires it, --no-sketch skips it. A --from-workdir sketch of the same genome is copied. |

### reads

Select sequencing runs from ENA/SRA by taxon or accession (writes reads.tsv).

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory (created). |
| `-tf`, `--target-family` |  | Restrict the selection to this family. Only the most specific of -tf/-tg/-ts given is used; they are not combined. |
| `-tg`, `--target-genus` |  | Restrict the selection to this genus. Only the most specific of -tf/-tg/-ts given is used; they are not combined. |
| `-ts`, `--target-species` |  | Restrict the selection to this species. Only the most specific of -tf/-tg/-ts given is used; they are not combined. |
| `--accession` |  | A run (SRR/ERR/DRR), sample (SAMN.., SRS..) or study (PRJNA.., SRP..) accession to include (repeatable). |
| `--accession-file` |  | File of accessions, one per line (# comments allowed). |
| `--platform` | `any` | Keep runs of one platform: any, illumina, ont or pacbio. |
| `--max-runs` |  | Keep at most N runs, the largest by bases. |
| `--min-bases` | `0` | Drop runs with fewer sequenced bases than this. |
| `--max-bases` |  | Drop runs with more sequenced bases than this (a guard against whole-host libraries, which would assemble into a host-dominated genome). |
| `--drop-selection` | `MDA` | Drop runs whose ENA library selection is this value (repeatable; default MDA, whole-genome amplification). Pass 'none' to keep every selection. |
| `--one-per-sample`, `--all-runs` | on | Keep the best run of each sample (a long-read run with enough bases, else the largest run), or every run. |

### assemble

Fetch and assemble the selected runs; write genomes/ and selection.tsv.

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory. |
| `--assembler` | `auto` | Assembler: auto, flye, shovill, skesa. |
| `-t`, `--threads` | `16` | Threads for the external tool (default 16, or the CPU limit when lower). |
| `--jobs` |  | Runs assembled concurrently, threads split across them (default 2, or 1 when a long-read run is pending). |
| `--memory-gb` | `16` | Memory hint per assembly, in GB, for tools that cap RAM; also the budget that caps concurrent classifier gathers (about 0.6 GB each). |
| `--min-contig-length` | `500` | Drop contigs shorter than this many bases. |
| `--polisher` | `auto` | Polisher for long-read assemblies: none, auto, medaka, racon. |
| `--polish-rounds` | `1` | Polishing rounds (racon; medaka runs one). |
| `--outgroup` |  | A FASTA file to set aside as the outgroup for rooting. |
| `--max-runs` |  | Keep at most N runs with a finished assembly, in reads.tsv order (largest first); finished runs count toward N and are kept. The other runs are listed as deferred in excused_runs.tsv; a later call with a larger N assembles them. |
| `--append` | off | Add the assemblies to a working directory that already holds a selection (metadata and genome, or ingest) instead of replacing it. |
| `--keep-reads` | off | Keep the downloaded FASTQ files after assembling. |
| `--keep-files` | off | Keep each run's assembler scratch directory. |
| `--checkm2-db` |  | CheckM2 DIAMOND database; enables quality scoring (or set CHECKM2DB). |
| `--min-completeness` | `50.0` | CheckM2 completeness floor. |
| `--max-contamination` | `10.0` | CheckM2 contamination ceiling. |
| `--classifier` | `auto` | Classifier: none, auto, sourmash. |
| `--gtdb-sketch` |  | GTDB sourmash sketch database (.sig.zip); enables classification (or set REPGENR_GTDB_SKETCH). |
| `--gtdb-lineages` |  | The lineages CSV published with the sketch (or set REPGENR_GTDB_LINEAGES). |
| `--tool-arg` |  | Assembler tuning as key=value (repeatable), e.g. mode=nano-raw. |
| `--sketch`, `--no-sketch` | auto | Write a sourmash sketch of each genome to sketches/ (k=21,31,51, scaled=1000). Default: when sourmash can run; --sketch requires it, --no-sketch skips it. With --append only the new genomes are sketched. |
| `--reads-sketch`, `--no-reads-sketch` | auto | Sketch each run's reads with sourmash (k=21,31,51, scaled=1000, with abundances) into assemblies/<run>/reads.sig.zip while the run is assembled; never in sketches/. Default: when sourmash can run; --reads-sketch requires it, --no-reads-sketch skips it. |
| `--screen-reads` | off | Screen each run's reads sketch before assembling it: a run is excused when its GTDB genus differs from the submitted genus (taxon_mismatch), when the top species holds less than --screen-min-fraction of the reads (low_match_fraction), or when an accepted run of the same biosample contains it (duplicate_isolate). Needs sourmash, --gtdb-sketch and --gtdb-lineages. |
| `--screen-min-fraction` | `0.3` | With --screen-reads: minimum abundance-weighted fraction of the reads in the top GTDB species. |
| `--screen-dup-ani` | `0.999` | With --screen-reads: ANI estimate from which a run contained in an earlier accepted run of the same biosample is excused as a duplicate. |

## Core stages

### dereplicate

Cluster genomes by ANI and select representatives.

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory. |
| `--tool` | `skder` | auto, drep, galah, skder, sourmash. |
| `-pani`, `--primary-ani` | `0.9` | Primary (pre-clustering) ANI threshold in (0, 1]. |
| `-sani`, `--secondary-ani` | `0.99` | Secondary (final cluster) ANI threshold in (0, 1]. |
| `-af`, `--aligned-fraction` | `0.5` | Minimum aligned fraction in (0, 1] for a pair to be compared. |
| `-t`, `--threads` | `16` | Threads for the external tool (default 16, or the CPU limit when lower). |
| `-s`, `--process-size` |  | Chunk size; when set and exceeded, two-stage chunking runs for any tool. |
| `-p`, `--num-processes` | `0` | Parallel stage-1 chunk workers (threads split across them). 0 = auto (~threads/4, capped by cores). |
| `--pre-primary-ani` |  | Stage-1 (intra-chunk) primary ANI; defaults to --primary-ani. |
| `--pre-secondary-ani` |  | Stage-1 (intra-chunk) secondary ANI; defaults to --secondary-ani. |
| `--reduce` | `none` | Taxonomy-aware reduction after ANI: none, species, or genus (one representative per taxon). |
| `--target-reps` | `0` | Target representative count: search --secondary-ani to land near it (0 = off; re-runs dereplication per search step). |
| `--virus` | off | Pass virus-tuned parameters to dRep (--tool drep); the other tools do not read it. |
| `--tool-arg` |  | Tool tuning as key=value (repeatable), e.g. mode=greedy. |
| `--allow-incomplete` | off | Proceed with a warning when the input genome set is incomplete. |
| `--keeper` | `quality` | Representative choice per cluster: quality (manifest CheckM values and genome N50), gtdb (a GTDB species representative first, then quality) or tool (adapter's own). |
| `--refine` |  | Seed of a lineage to keep at finer resolution (filename, stem or accession; repeatable): every cluster holding a seed is dereplicated again at --refine-ani. |
| `--refine-ani` |  | Secondary ANI for the clusters of the --refine seeds; above --secondary-ani. |

### derep-refine

Keep more genomes inside a lineage: dereplicate the seeds' clusters again at a finer ANI.

Refines the dereplication on disk with the tool and settings of its
record, without repeating it. Store the current result first with
derep-stock if it should be kept.

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory. |
| `--refine` | required | Seed of a lineage to keep at finer resolution (filename, stem or accession; repeatable): every cluster holding a seed is dereplicated again at --refine-ani. |
| `--refine-ani` | required | Secondary ANI for the clusters of the --refine seeds; above --secondary-ani. |
| `-t`, `--threads` | `16` | Threads for the external tool (default 16, or the CPU limit when lower). |

### annotate

Annotate the genome set once (GFF3 per genome, cached) for the core-gene route.

On this command --tool selects the annotator.

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory. |
| `--tool` | `bakta` | Annotator: bakta (default; needs its database, --bakta-db or BAKTA_DB) or prokka. |
| `--all-genomes` | off | Use all genomes, not only the representatives. |
| `--bakta-db` |  | bakta database directory (bakta_db download --output DIR --type light\|full); BAKTA_DB when unset. |
| `-t`, `--threads` | `16` | Threads for the external tool (default 16, or the CPU limit when lower). |
| `--tool-arg` |  | Annotator tuning as key=value (repeatable). |
| `--allow-incomplete` | off | Proceed with a warning when the input genome set is incomplete. |
| `-p`, `--num-processes` | `0` | Genomes annotated at once (threads split across them). 0 = auto (~threads/4, capped by cores). |

### snptype

Call SNPs and build a core-SNP alignment.

On this command --tool selects the SNP typer, whereas on the dereplication
commands it selects the dereplicator.

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory. |
| `--tool` | `simple` | SNP typer: parsnp, simple, ska2, snippy. |
| `--reference` |  | Reference genome: a filename, 'medoid' (the default: the genome nearest to all others by sketch ANI, when the genomes are sketched and sourmash is at hand) or 'first'. |
| `--all-genomes` | off | Use all genomes, not only the representatives. |
| `--mask` | `none` | Recombination masking of the SNP alignment: none, gubbins. |
| `--core-fraction` |  | Keep only the variable columns where at least this fraction of the genomes have a called base (a soft core, 0 to 1). Unset keeps every variable column and the log reports what 0.95 would keep. |
| `-t`, `--threads` | `16` | Threads for the external tool (default 16, or the CPU limit when lower). |
| `--tool-arg` |  | Tool tuning as key=value (repeatable). |
| `--allow-incomplete` | off | Proceed with a warning when the input genome set is incomplete. |

### phylo

Build a phylogenetic tree from an alignment, SNP alignment, or directly.

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory. |
| `--treebuilder` | `iqtree` | auto, fasttree, iqtree, mashtree, raxmlng, ska2, skani, sourmash, veryfasttree. |
| `--msa-source` | `aligner` | aligner, snptype or coregene (a concatenated core-gene alignment from the annotations the annotate stage wrote). |
| `--aligner` | `progressivemauve` | cactus, progressivemauve, sibeliaz. |
| `--snptyper` | `simple` | SNP typer: parsnp, simple, ska2, snippy. |
| `--all-genomes` | off | Use all genomes, not only the representatives. |
| `--no-outgroup` | off | Do not root with an outgroup. |
| `--root` |  | How phylo roots tree/tree.nwk: outgroup (on the branch to the outgroup), midpoint, or none (as the tree builder wrote it). Default: outgroup when an outgroup is staged and --no-outgroup is not given, else none. |
| `-B`, `--bootstrap` | `0` | Bootstrap replicates. 0 turns bootstrapping off; IQ-TREE needs at least 1000 when it is on. |
| `--reference` |  | Reference genome: a filename, 'medoid' (the default: the genome nearest to all others by sketch ANI, when the genomes are sketched and sourmash is at hand) or 'first'. |
| `--aligner-arg` |  | Aligner, SNP typer, masker or tree builder tuning as key=value (repeatable), e.g. kmer=15 (sibeliaz), seed_weight=11 (progressivemauve) or model=GTR+G (iqtree, raxmlng). |
| `-t`, `--threads` | `16` | Threads for the external tool (default 16, or the CPU limit when lower). |
| `--mask` | `none` | Recombination masking of the SNP alignment: none, gubbins. Needs --msa-source snptype. |
| `--core-fraction` |  | Keep only the variable columns where at least this fraction of the genomes have a called base (a soft core, 0 to 1). Unset keeps every variable column and the log reports what 0.95 would keep. |
| `--coregene` | `panaroo` | Core-gene alignment builder for --msa-source coregene: panaroo, ppanggolin, scarap or orthofinder. Tuning goes through --aligner-arg. panaroo: threshold= is the identity at which it groups genes (0.98; lower it, 0.95 then 0.90, when the core comes out small on a divergent set), and also filtered=false, merge_paralogs, remove_invalid_genes, family_threshold, len_dif_percent. ppanggolin: identity= and coverage= of its protein clustering (0.8 each), partition=softcore\|core\|persistent (softcore, at --core-threshold), single_copy=false, dup_margin, mode, nb_of_partitions, translation_table. scarap: method (FH), seeds (100), max_core_genes (0, no cap), core_prefilter (the smaller of 0.90 and --core-threshold). orthofinder: search=diamond\|diamond_ultra_sens\|blastp\|mmseqs (diamond), inflation (MCL, 1.2). |
| `--core-threshold` |  | With --msa-source coregene: a gene is core when present in at least this fraction of the genomes (0.95). |
| `--clean-mode` |  | With --msa-source coregene: Panaroo's cleaning of annotation errors, strict (default), moderate or sensitive. |
| `--allow-incomplete` | off | Proceed with a warning when the input genome set is incomplete. |

### tree2tax

Emit FlexTaxD-compatible taxonomy relations from the tree.

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory. |
| `--node-basename` |  | Name internal nodes <basename><n>. Without it, internal nodes receive names derived from a hash of their descendant leaves. |
| `--naming` | `leafset` | How internal nodes are named without --node-basename: 'leafset' (a hash of the whole leaf set) or 'lca' (a hash of the smallest leaf label under each child, so most names survive genomes added with later-sorting accessions). |
| `--root-name` | `root` | Label of the top node. FlexTaxD expects 'root'. |
| `--remove-outgroup` | off | Leave the outgroup out of the taxonomy after rooting. |
| `--include-dereplicated`, `--no-include-dereplicated` | on | List redundant genomes under their representative in the taxonomy. |
| `--collapse-support` |  | Merge nodes whose support is below this fraction into their parent. |
| `--collapse-length` |  | Merge nodes whose branch is shorter than this length into their parent. |
| `--collapse-diameter` |  | Merge every clade whose leaf-to-leaf distance is within this value into one node, keeping the fewest such clades (in the tree's own units). |
| `--ranks` | off | Add a rank column to tree2tax.tsv: 'strain' on every genome, 'no rank' on internal nodes unless --rank-red names them. |
| `--rank-red` |  | With --ranks: NAME=RED bands, e.g. species=0.95,genus=0.8; an internal node takes the rank with the largest threshold not above its relative evolutionary divergence, once per lineage. Needs branch lengths and a rooted tree. |

## Inspect a dereplication

### census

Count the genera, species and samples under a taxon (read-only).

Without -wd, queries a GTDB family (-tf) or genus (-tg) through the GTDB
API or a GTDB metadata table, optionally with the ENA sequencing runs
(--runs), or a viral taxon through NCBI Virus (--viral --target; metadata
only, no sequences). With -wd, counts the candidates the entry stage found
or the genomes it selected, by source; -tf and -tg narrow the count. A
family is counted per genus and a genus per species. No stage is recorded
and nothing is written to a working directory. Exits 2 for a missing taxon
or an unsupported combination, 3 when a request fails or -wd does not exist.

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` |  | Count the candidates (after metadata or vmetadata) or the selection (after genome, vgenome, ingest or assemble) of this working directory; nothing is written to it. |
| `-tf`, `--target-family` |  | Count this family: one row per genus. |
| `-tg`, `--target-genus` |  | Count this genus: one row per species. |
| `--viral` | off | Count a viral taxon (NCBI Virus) instead of a GTDB taxon. |
| `--target` |  | With --viral: the virus taxon, as vmetadata takes it (e.g. picornaviridae). |
| `--source` |  | api (GTDB API, default) or table (GTDB metadata table, also 'tsv'); with --viral, ncbi_virus (default). bvbrc is counted only from a vmetadata workdir. |
| `-r`, `--release` |  | GTDB release of the table source, e.g. 232.0. |
| `--gtdb-version` |  | GTDB table of the table source: bac120 (default) or ar53. |
| `--metadata-path` |  | Read this GTDB metadata table (.tsv.gz) instead of downloading the release table; with -wd, the table to count the candidates from. |
| `--runs` | off | Also count the ENA whole-genome sequencing runs under the taxon: runs, biosamples and runs per platform (bacterial taxa only). |
| `--host` |  | With --viral: only records from this host species. |
| `--complete-only` | off | With --viral: only sequences marked complete. |
| `--released-after` |  | With --viral: only records released after this date (MM/DD/YYYY). |
| `--tsv` |  | Also write the rows to this TSV file (format in docs/output.md). |
| `--json` | off | Print one JSON object (taxon, mode, totals, rows) on stdout instead of the table. |

### glance

Quick all-vs-all ANI overview (dendrogram + similarity plots).

Needs dRep or sourmash (on the PATH, or via the container backend) and
no dereplication; compares every genome in genomes/ and writes
glance_clustering_dendrogram.pdf and two ANI plots,
glance_MASH_ANI_similarity_boxplot.png and
glance_MASH_ANI_similarity_histogram.png. dRep reports Mash ANI;
sourmash reports the k-mer ANI estimate that dereplicate --tool
sourmash thresholds.

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory. |
| `--tool` | `auto` | Dereplicator with comparison support: auto, drep, sourmash. auto uses dRep when it can run (on the PATH or via the container backend), otherwise sourmash. |
| `-t`, `--threads` | `16` | Threads for the external tool (default 16, or the CPU limit when lower). |
| `--plot-max` | `1.0` | Upper bound of the ANI values plotted, as a fraction from 0 to 1. |
| `--plot-min` | `0.0` | Lower bound of the ANI values plotted, as a fraction from 0 to 1. |
| `--keep-files` | off | Keep the comparison tool's working directory glance_wd/. |

### cluster-summary

Regenerate derep/cluster_summary.tsv (size, species, keeper quality per cluster).

n_members counts the genomes under the representative and excludes it;
n_genomes includes it. Species come from the manifest taxonomy, or from
canonical filenames when the manifest has none, and at most five are
listed.

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory. |

### derep-unpack

Explode clusters into one directory per representative.

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory. |
| `--no-representant` | off | Leave the representative out of its cluster directory. |

### derep-stock

Store, load, list or delete named dereplication runs.

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory. |
| `--action` | required | list, pack, unpack or delete. |
| `--name` |  | Run name for pack/unpack/delete: up to 100 letters, digits, '.', '_' or '-', starting with a letter or digit. |

### sketch

Write the sourmash sketch of each genome to sketches/.

One file per genome, sketches/<name>.sig.zip, with DNA signatures at
k=21, 31 and 51 (scaled=1000) named after the genome. Only missing and
stale sketches are written (stale: the FASTA or the parameters changed);
sketches of genomes no longer in the set are removed. repgenr --force
sketch writes every sketch again. Needs sourmash (on the PATH, or via the
container backend).

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory. |
| `-t`, `--threads` | `16` | Threads for the external tool (default 16, or the CPU limit when lower). |

## Environment and diagnostics

### list-tools

List the available pluggable tools in each family.

A tool that declares a recommended scale is shown as 'name (up to N
genomes)'; auto-selection and the scale warnings use the same limit.
The last line names the dereplicators that glance can run.
With --check, every adapter's required binaries are looked up (version
floors included) and reported per tool, so an environment can be
verified before a run without a working directory. Under a container
backend each line names where the tool runs: '[image <ref>]' or
'[host]'; --images adds, for each tool that passed, whether its images
(secondary ones such as racon's minimap2 included) are present
locally. A version query that does not answer within 8 s is stopped,
and the version is shown as unknown unless the tool's conda package
record names one. --check alone
always exits 0, since a host that has only some families installed is
normal; --check --strict exits 4 when any adapter is missing or errored,
or 5 when any plugin failed to load, so a script can verify an
environment.

| option | default | description |
|---|---|---|
| `--check` | off | Run each adapter's preflight and report whether its binaries are present. |
| `--strict` | off | With --check, exit 4 when an adapter is missing or errored and 5 when a plugin is broken, after the full listing. |
| `--images` | off | With --check under a container backend, also report whether each tool's image is present locally (docker image inspect or the Singularity .sif cache; nothing is pulled). |

### doctor

Verify a workdir's outputs against its records (read-only health check).

`status` lists each stage as done, stale or interrupted; `doctor` also
checks the outputs themselves: missing, corrupt or untracked genomes,
dangling links, manifest drift, representative/cluster mismatches,
truncated tree and tree2tax tables, missing deliverables, stages whose
inputs changed since completion, and leftover temp files.
Exits 0 when only warnings are found (a stale stage re-runs on its next
invocation), 7 when any failure is found (including a malformed
repgenr.yaml or a check that could not complete), 3 when the workdir
does not exist, and 1 only on an unexpected error.

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory. |
| `--quick` | off | Skip reading the first bytes of each genome file (the FASTA content check), the slowest check on large genome sets. Links, missing and untracked genomes are still checked. |
| `--json` | off | Print one versioned JSON object on stdout instead of the text report (schema in docs/output.md); exit codes are unchanged. |

### versions

Print the external-tool versions recorded in a workdir's repgenr.yaml.

One 'tool: version' line per tool. A tool that stages recorded with
different versions (an image in one, the host binary in another) is listed
once per stage as 'tool (stage): version'. A containerized tool's version
is its image reference. A stage that did not finish is named on stderr.

Lets the Nextflow bridge modules (which run a full stage in a scratch
workdir) surface the resolved tool versions into versions.yml.

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory. |
| `--versions-out` |  | Write a versions.yml fragment here instead of stdout. |

## Nextflow data-channel steps

### genome-fetch

Download genomes listed in a selection.tsv (stateless data-channel step).

| option | default | description |
|---|---|---|
| `--selection` | required | selection.tsv from the metadata stage. |
| `-o`, `--out` | required | Output dir for downloaded genomes. |
| `--keep-files` | off | Keep download and scratch intermediates. |
| `--versions-out` |  | Write resolved tool versions (YAML fragment) here. |

### dereplicate-chunk

Dereplicate one chunk of genomes (scatter step; writes a chunk result dir).

| option | default | description |
|---|---|---|
| `--genomes-fofn` | required | File of genome FASTA paths, one per line. |
| `-o`, `--out` | required | Output directory for the chunk result. |
| `--tool` | `skder` | drep, galah, skder, sourmash. |
| `-pani`, `--primary-ani` | `0.9` | Primary (pre-clustering) ANI threshold in (0, 1]. |
| `-sani`, `--secondary-ani` | `0.99` | Secondary (final cluster) ANI threshold in (0, 1]. |
| `-af`, `--aligned-fraction` | `0.5` | Minimum aligned fraction in (0, 1] for a pair to be compared. |
| `-t`, `--threads` | `16` | Threads for the external tool (default 16, or the CPU limit when lower). |
| `--virus` | off | Pass virus-tuned parameters to dRep (--tool drep); the other tools do not read it. |
| `--tool-arg` |  | Tool tuning as key=value (repeatable), e.g. mode=greedy. |
| `--selection-tsv` |  | selection.tsv with quality columns; enables quality-aware representatives. |
| `--keeper` | `quality` | Representative choice when --selection-tsv is given: quality (selection.tsv CheckM values and genome N50), gtdb (GTDB species representative first, then quality) or tool (adapter's own pick). |
| `--sketches-dir` |  | A sketches/ directory of <record name>.sig.zip files (k=21, 31, 51; scaled=1000) made from these genomes. A sourmash tool reads them instead of sketching; other tools ignore it. |
| `--versions-out` |  | Write resolved tool versions (YAML fragment) here. |

### dereplicate-merge

Dereplicate the union of chunk representatives (gather step).

| option | default | description |
|---|---|---|
| `-o`, `--out` | required | Output dir for the merged result. |
| `--chunk-dir` |  | A chunk result directory (repeatable). |
| `--chunk-fofn` |  | File listing chunk result directories, one per line. |
| `--tool` | `skder` | drep, galah, skder, sourmash. |
| `-pani`, `--primary-ani` | `0.9` | Primary (pre-clustering) ANI threshold in (0, 1]. |
| `-sani`, `--secondary-ani` | `0.99` | Secondary (final cluster) ANI threshold in (0, 1]. |
| `-af`, `--aligned-fraction` | `0.5` | Minimum aligned fraction in (0, 1] for a pair to be compared. |
| `-t`, `--threads` | `16` | Threads for the external tool (default 16, or the CPU limit when lower). |
| `--virus` | off | Pass virus-tuned parameters to dRep (--tool drep); the other tools do not read it. |
| `--tool-arg` |  | Tool tuning as key=value (repeatable), e.g. mode=greedy. |
| `--selection-tsv` |  | selection.tsv with quality columns; enables quality-aware representatives. |
| `--keeper` | `quality` | Representative choice when --selection-tsv is given: quality (selection.tsv CheckM values and genome N50), gtdb (GTDB species representative first, then quality) or tool (adapter's own pick). |
| `--reduce` | `none` | Taxonomy-aware reduction after the merge: none, species, or genus (one representative per taxon; taxonomy from --selection-tsv or the filenames). |
| `--target-reps` | `0` | Target representative count: search --secondary-ani of the merge pass to land near it (0 = off; re-runs the merge per search step). |
| `--sketches-dir` |  | A sketches/ directory of <record name>.sig.zip files (k=21, 31, 51; scaled=1000) made from these genomes. A sourmash tool reads them instead of sketching; other tools ignore it. |
| `--versions-out` |  | Write resolved tool versions (YAML fragment) here. |

### annotate-build

Annotate one genome (stateless data-channel step; writes <out>/<name>/<name>.gff3).

| option | default | description |
|---|---|---|
| `--genome` | required | Genome FASTA to annotate (gz accepted). |
| `-o`, `--out` | required | Output directory; the record's files go under <out>/<name>/. |
| `--record-name` |  | Record name (prefix of the outputs); the genome's by default. |
| `--tool` | `bakta` | Annotator: bakta (default; needs its database, --bakta-db or BAKTA_DB) or prokka. |
| `--bakta-db` |  | bakta database directory (bakta_db download --output DIR --type light\|full); BAKTA_DB when unset. |
| `-t`, `--threads` | `16` | Threads for the external tool (default 16, or the CPU limit when lower). |
| `--tool-arg` |  | Annotator tuning as key=value. |
| `--versions-out` |  | Write resolved tool versions (YAML fragment) here. |

### annotate-index

Index gathered annotate-build outputs (stateless data-channel step; writes
<dir>/annotations.tsv, which phylo-build --annotations-dir reads).

| option | default | description |
|---|---|---|
| `--annotations-dir` | required | Directory holding the record directories annotate-build wrote (<dir>/<name>/<name>.gff3). |
| `--tool` | `bakta` | Annotator assumed for a record without annotation.json (output of an earlier version): bakta or prokka. Each other record is indexed with the annotator it names. |

### phylo-build

Build a phylogeny from a genomes directory (stateless data-channel step).

Here --outgroup-accession takes a file that names the accession, not the accession itself.

| option | default | description |
|---|---|---|
| `--genomes-dir` | required | Directory of genome FASTA files to build the tree from. |
| `-o`, `--out` | required | Output dir (writes tree/tree.nwk). |
| `--outgroup-dir` |  | Directory holding the outgroup genome file(s). |
| `--outgroup-accession` |  | File naming the outgroup accession. |
| `--treebuilder` | `iqtree` | auto, fasttree, iqtree, mashtree, raxmlng, ska2, skani, sourmash, veryfasttree. |
| `--msa-source` | `aligner` | aligner, snptype or coregene (a concatenated core-gene alignment from the annotations the annotate stage wrote). |
| `--aligner` | `progressivemauve` | cactus, progressivemauve, sibeliaz. |
| `--snptyper` | `simple` | SNP typer: parsnp, simple, ska2, snippy. |
| `--no-outgroup` | off | Do not root with an outgroup. |
| `--root` |  | How phylo roots tree/tree.nwk: outgroup (on the branch to the outgroup), midpoint, or none (as the tree builder wrote it). Default: outgroup when an outgroup is staged and --no-outgroup is not given, else none. |
| `-B`, `--bootstrap` | `0` | Bootstrap replicates. 0 turns bootstrapping off; IQ-TREE needs at least 1000 when it is on. |
| `--core-fraction` |  | Keep only the variable columns where at least this fraction of the genomes have a called base (a soft core, 0 to 1). Unset keeps every variable column and the log reports what 0.95 would keep. |
| `--coregene` | `panaroo` | Core-gene alignment builder for --msa-source coregene: panaroo, ppanggolin, scarap or orthofinder. Tuning goes through --aligner-arg. panaroo: threshold= is the identity at which it groups genes (0.98; lower it, 0.95 then 0.90, when the core comes out small on a divergent set), and also filtered=false, merge_paralogs, remove_invalid_genes, family_threshold, len_dif_percent. ppanggolin: identity= and coverage= of its protein clustering (0.8 each), partition=softcore\|core\|persistent (softcore, at --core-threshold), single_copy=false, dup_margin, mode, nb_of_partitions, translation_table. scarap: method (FH), seeds (100), max_core_genes (0, no cap), core_prefilter (the smaller of 0.90 and --core-threshold). orthofinder: search=diamond\|diamond_ultra_sens\|blastp\|mmseqs (diamond), inflation (MCL, 1.2). |
| `--core-threshold` |  | With --msa-source coregene: a gene is core when present in at least this fraction of the genomes (0.95). |
| `--clean-mode` |  | With --msa-source coregene: Panaroo's cleaning of annotation errors, strict (default), moderate or sensitive. |
| `--annotations-dir` |  | With --msa-source coregene: an annotations/ directory the annotate stage or annotate-build wrote (its index names the records); without it the step annotates the genomes itself with --annotator. |
| `--annotator` | `bakta` | Annotator: bakta (default; needs its database, --bakta-db or BAKTA_DB) or prokka. |
| `--bakta-db` |  | bakta database directory (bakta_db download --output DIR --type light\|full); BAKTA_DB when unset. |
| `--reference` |  | Reference genome: a filename, 'medoid' (the default: the genome nearest to all others by sketch ANI, when the genomes are sketched and sourmash is at hand) or 'first'. |
| `--aligner-arg` |  | Aligner, SNP typer, masker or tree builder tuning as key=value (repeatable), e.g. kmer=15 (sibeliaz), seed_weight=11 (progressivemauve) or model=GTR+G (iqtree, raxmlng). |
| `-t`, `--threads` | `16` | Threads for the external tool (default 16, or the CPU limit when lower). |
| `--mask` | `none` | Recombination masking of the SNP alignment: none, gubbins. Needs --msa-source snptype. |
| `--msa-only` | off | Build the alignment and stop, writing msa.fasta (for a separate tree step). |
| `--msa` |  | Build the tree from this alignment instead of constructing one. |
| `--sketches-dir` |  | A sketches/ directory of <record name>.sig.zip files (k=21, 31, 51; scaled=1000) made from these genomes. A sourmash tool reads them instead of sketching; other tools ignore it. |
| `--versions-out` |  | Write resolved tool versions (YAML fragment) here. |

### tree2tax-relations

Emit FlexTaxD relations from a tree (stateless data-channel step).

Here --outgroup-accession takes a file that names the accession, not the accession itself.

| option | default | description |
|---|---|---|
| `--tree` | required | Rooted/unrooted tree in Newick (tree.nwk). |
| `-o`, `--out` | required | Output dir (writes tree2tax.tsv + genomes_map.tsv). |
| `--clusters` |  | derep clusters.tsv (for --include-dereplicated). |
| `--outgroup-dir` |  | Directory holding the outgroup genome file(s). |
| `--outgroup-accession` |  | File naming the outgroup accession. |
| `--segments` |  | segments.tsv from vgenome --group-segments: member segment accessions are listed under their isolate's leaf. |
| `--node-basename` |  | Name internal nodes <basename><n>. Without it, internal nodes receive names derived from a hash of their descendant leaves. |
| `--naming` | `leafset` | How internal nodes are named without --node-basename: 'leafset' (a hash of the whole leaf set) or 'lca' (a hash of the smallest leaf label under each child, so most names survive genomes added with later-sorting accessions). |
| `--root-name` | `root` | Label of the top node. FlexTaxD expects 'root'. |
| `--remove-outgroup` | off | Leave the outgroup out of the taxonomy after rooting. |
| `--no-outgroup` | off | Do not root with an outgroup. |
| `--keep-root` | off | Keep the tree's own root (a tree from phylo --root midpoint) instead of rerooting on the outgroup; the outgroup is still read for --remove-outgroup. |
| `--include-dereplicated`, `--no-include-dereplicated` | on | List redundant genomes under their representative in the taxonomy. |
| `--versions-out` |  | Write resolved tool versions (YAML fragment) here. |
| `--collapse-support` |  | Merge nodes whose support is below this fraction into their parent. |
| `--collapse-length` |  | Merge nodes whose branch is shorter than this length into their parent. |
| `--collapse-diameter` |  | Merge every clade whose leaf-to-leaf distance is within this value into one node, keeping the fewest such clades (in the tree's own units). |
| `--ranks` | off | Add a rank column to tree2tax.tsv: 'strain' on every genome, 'no rank' on internal nodes unless --rank-red names them. |
| `--rank-red` |  | With --ranks: NAME=RED bands, e.g. species=0.95,genus=0.8; an internal node takes the rank with the largest threshold not above its relative evolutionary divergence, once per lineage. Needs branch lengths and a rooted tree. |

### assemble-run

Fetch and assemble one run of a reads.tsv (stateless data-channel step).

| option | default | description |
|---|---|---|
| `--reads-tsv` | required | reads.tsv from the reads stage. |
| `--run` | required | The run accession (a row of reads.tsv) to assemble. |
| `-o`, `--out` | required | Output dir: contigs.fasta and assembly.ok, or excused_runs.tsv. |
| `--assembler` | `auto` | Assembler: auto, flye, shovill, skesa. |
| `-t`, `--threads` | `16` | Threads for the external tool (default 16, or the CPU limit when lower). |
| `--memory-gb` | `16` | Memory hint for the assembly, in GB, for tools that cap RAM. |
| `--min-contig-length` | `500` | Drop contigs shorter than this many bases. |
| `--polisher` | `auto` | Polisher for long-read assemblies: none, auto, medaka, racon. |
| `--polish-rounds` | `1` | Polishing rounds (racon; medaka runs one). |
| `--keep-reads` | off | Keep the downloaded FASTQ files after assembling. |
| `--keep-files` | off | Keep the assembler scratch directory. |
| `--tool-arg` |  | Assembler tuning as key=value (repeatable), e.g. mode=nano-raw. |
| `--reads-sketch`, `--no-reads-sketch` | auto | Sketch each run's reads with sourmash (k=21,31,51, scaled=1000, with abundances) into assemblies/<run>/reads.sig.zip while the run is assembled; never in sketches/. Default: when sourmash can run; --reads-sketch requires it, --no-reads-sketch skips it. |
| `--screen-reads` | off | Screen each run's reads sketch before assembling it: a run is excused when its GTDB genus differs from the submitted genus (taxon_mismatch), when the top species holds less than --screen-min-fraction of the reads (low_match_fraction), or when an accepted run of the same biosample contains it (duplicate_isolate). Needs sourmash, --gtdb-sketch and --gtdb-lineages. The duplicate check needs the other runs and is not made by this step. |
| `--screen-min-fraction` | `0.3` | With --screen-reads: minimum abundance-weighted fraction of the reads in the top GTDB species. |
| `--gtdb-sketch` |  | GTDB sourmash sketch database for --screen-reads (or set REPGENR_GTDB_SKETCH). |
| `--gtdb-lineages` |  | The lineages CSV published with the sketch (or set REPGENR_GTDB_LINEAGES). |
| `--versions-out` |  | Write resolved tool versions (YAML fragment) here. |

### genome-qc

Score (CheckM2) and classify a batch of assemblies (stateless data-channel step).

| option | default | description |
|---|---|---|
| `--assemblies` | required | Directory of assemble-run output dirs, one per run. |
| `-o`, `--out` | required | Output dir for quality.tsv and classification.tsv. |
| `-t`, `--threads` | `16` | Threads for the external tool (default 16, or the CPU limit when lower). |
| `--memory-gb` | `16` | Memory budget in GB; caps concurrent classifier gathers at about 0.6 GB each. |
| `--checkm2-db` |  | CheckM2 DIAMOND database; enables quality scoring (or set CHECKM2DB). |
| `--classifier` | `auto` | Classifier: none, auto, sourmash. |
| `--gtdb-sketch` |  | GTDB sourmash sketch database (.sig.zip); enables classification (or set REPGENR_GTDB_SKETCH). |
| `--gtdb-lineages` |  | The lineages CSV published with the sketch (or set REPGENR_GTDB_LINEAGES). |
| `--tool-arg` |  | Classifier tuning as key=value (repeatable). |
| `--versions-out` |  | Write resolved tool versions (YAML fragment) here. |

### reads-gather

Write the genome contract from per-run assemblies (stateless data-channel step).

| option | default | description |
|---|---|---|
| `--reads-tsv` | required | reads.tsv from the reads stage. |
| `--assemblies` | required | Directory of assemble-run output dirs, one per run. |
| `-o`, `--out` | required | Output dir for genomes/, selection.tsv and the stats tables. |
| `--qc` |  | genome-qc output dir (quality.tsv, classification.tsv), if it ran. |
| `--min-completeness` | `50.0` | CheckM2 completeness floor. |
| `--max-contamination` | `10.0` | CheckM2 contamination ceiling. |
