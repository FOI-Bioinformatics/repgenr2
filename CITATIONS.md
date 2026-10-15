# RepGenR2: Citations

## Pipeline

RepGenR2 builds representative-genome repositories and the associated FlexTaxD
taxonomy. If you use it, please cite this repository:
<https://github.com/FOI-Bioinformatics/repgenr2>.

## Workflow managers

- [Nextflow](https://doi.org/10.1038/nbt.3820)

  > Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C.
  > Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017.

## Tools

The stages below shell out to external tools; cite the ones your run uses.

- **Genome selection / download**
  - [GTDB](https://doi.org/10.1093/nar/gkab776) (bacterial taxonomy / metadata)
  - [BV-BRC](https://doi.org/10.1093/nar/gkac1003) (viral metadata)
  - [NCBI Datasets](https://doi.org/10.1093/nar/gkad924) (genome download)

- **Dereplication**
  - [skDER / skani](https://doi.org/10.1186/s13059-024-03184-z)
  - [Galah](https://doi.org/10.1186/s40168-021-01213-8)
  - [sourmash](https://doi.org/10.21105/joss.00027) and the
    [branchwater plugin](https://doi.org/10.21105/joss.06830)
  - [dRep](https://doi.org/10.1038/ismej.2017.126), which runs
    [CheckM](https://doi.org/10.1101/gr.186072.114) for genome quality,
    [Prodigal](https://doi.org/10.1186/1471-2105-11-119) for gene calling,
    [Mash](https://doi.org/10.1186/s13059-016-0997-x) for its primary clustering
    and [FastANI](https://doi.org/10.1038/s41467-018-07641-9) or skani for the
    secondary comparisons
  - [Mash](https://doi.org/10.1186/s13059-016-0997-x) also underlies Mashtree
    and the `glance` comparison

- **Alignment and phylogeny**
  - [progressiveMauve](https://doi.org/10.1371/journal.pone.0011147)
  - [SibeliaZ](https://doi.org/10.1038/s41467-020-19777-8)
  - [Cactus / Minigraph-Cactus](https://doi.org/10.1038/s41587-023-01793-w) and
    the [HAL toolkit](https://doi.org/10.1093/bioinformatics/btt128) (`hal2maf`)
  - [Mashtree](https://doi.org/10.21105/joss.01762)
  - [RAxML-NG](https://doi.org/10.1093/bioinformatics/btz305)
  - [IQ-TREE 2](https://doi.org/10.1093/molbev/msaa015)
  - [FastTree 2](https://doi.org/10.1371/journal.pone.0009490)

- **SNP typing**
  - [minimap2](https://doi.org/10.1093/bioinformatics/bty191),
    [SAMtools / BCFtools](https://doi.org/10.1093/gigascience/giab008)
  - [Snippy](https://github.com/tseemann/snippy)
  - [ParSNP / Harvest](https://doi.org/10.1186/s13059-014-0524-x)
  - [SKA2](https://doi.org/10.1101/gr.279449.124) (split k-mer analysis; Derelle et al. 2024, Genome Research)
  - [Gubbins](https://doi.org/10.1093/nar/gku1196) (recombination masking)

- **Annotation (core-gene route)**
  - [Bakta](https://doi.org/10.1099/mgen.0.000685) (Schwengers et al. 2021, Microbial Genomics)
  - [Prokka](https://doi.org/10.1093/bioinformatics/btu153) (Seemann 2014, Bioinformatics)
  - [Panaroo](https://doi.org/10.1186/s13059-020-02090-4) (Tonkin-Hill et al. 2020, Genome Biology)
  - [MAFFT](https://doi.org/10.1093/molbev/mst010) (Katoh and Standley 2013, Molecular Biology and Evolution)

- **Read assembly (reads chain)**
  - [ENA](https://doi.org/10.1093/nar/gkae975) (sequencing-run discovery and FASTQ files)
  - [SKESA](https://doi.org/10.1186/s13059-018-1540-z)
  - [Shovill](https://github.com/tseemann/shovill), which runs
    [SPAdes](https://doi.org/10.1089/cmb.2012.0021)
  - [Flye](https://doi.org/10.1038/s41587-019-0072-8)
  - [medaka](https://github.com/nanoporetech/medaka)
  - [racon](https://doi.org/10.1101/gr.214270.116)
  - [CheckM2](https://doi.org/10.1038/s41592-023-01940-w) (assembly quality)
  - [sourmash gather](https://doi.org/10.1101/2022.01.11.475838) against the
    [GTDB](https://doi.org/10.1093/nar/gkab776) sketch database (classification)

- **Taxonomy**
  - [FlexTaxD](https://doi.org/10.1093/bioinformatics/btab372)

## Software packaging

- [BioContainers](https://doi.org/10.1093/bioinformatics/btx192)
- [Seqera Wave](https://seqera.io/wave/)
- [Anaconda / conda-forge / Bioconda](https://doi.org/10.1038/s41592-018-0046-7)
