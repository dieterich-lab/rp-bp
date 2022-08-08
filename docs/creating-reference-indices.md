# More about creating reference indices

The [usage instructions](usage-instructions.html) give the necessary commands to
create the reference indices required for Rp-Bp. This document describes how to
download or otherwise find the necessary files and sequences.

<a id='toc'></a>

- [Reference genome sequence](#ref-genome)
- [Reference annotations](#ref-annotations)
- [The "ribosomal" index](#ribosomal-index)
- [Using _de novo_ assembly](#de-novo-assembly)

---

<a name="ref-genome"></a>

## Reference genome sequence

While the most recent version available of the reference genome sequence should
be used, we do not recommend to include things like haplotypes in the reference
for Rp-Bp. First, Rp-Bp does not treat these differently than "normal"
chromosomes. Second, the inclusion of these can greatly increase the alignment
step (with STAR). Most importantly, though, if these extra sequences
substantially overlap the main reference, reads which align to both will be
called as multimappers and discarded.

As mentioned in the [usage instructions](usage-instructions.html), the "primary
assembly" file from Ensembl contains the appropriate sequences and identifiers.
The "top level" Ensembl genome assembly includes haplotype information. Please
see [Ensembl](http://www.ensembl.org/info/genome/genebuild/assembly.html) for
more information about the differences between assemblies.

<a name="ref-annotations"></a>

## Reference annotations

Following standard conventions, the genomic locations are taken to be base-1 and
inclusive. Again, we recommend using the most recent version of the annotations, and in all cases the annotations **must match** the version of the reference genome sequence. At least the
`exon` features must be present, and the transcript identifiers (`transcript_id`
attribute for GTF) must match for exons from the same transcript. The `CDS` features
must also be present in the annotation file. Other features such as `gene`, `start_codon` and `stop_codon`
may be present in the file, but they are **not** used for extracting ORFs.

### A note regarding the annotated coding sequence (CDS)

There seem to be no clear consensus if the `STOP codon` should be part of the coding region.
The **GTF2 (GTF)** format specifications suggest to exclude it, and this is the default in Rp-Bp
when the annotation file as a .gtf extension. _Rp-Bp's ORF labeling algorithm
assumes that annotated CDSs include the START codon but exclude the STOP
codon_. If the annotations follow the [GFF3 format specifications](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md), then the STOP codon _is included_ and will be automatically removed during processing.
If unsure, this can be verified by loading the annotations and sequence into a
genome viewer such as IGV. If there is a problem with the annotations, it is likely to manifest
as many `within` and `suspect` ORFs while very few `canonical` ORFs predicted as translated.

**N.B** If the annotation file has a .gff extension, it will be treated according to the GFF3 format specifications. _Note that there is currently no "check" based on the format, and the selection is made solely based on the extension, as read from the configuration file._

**N.B.** In case CDS annotations are _not_ available, Rp-Bp can still be run as
normal. The only difference is that the ORF labels will not be meaningful!

<a name="ribosomal-index"></a>

## The "ribosomal" index

The ribosomal filtering step is not as imperative for Rp-Bp as for some other
types of analysis (such as anything based on RPKM); this is because reads which
align to these types of sequences tend to be multimappers and are discarded
anyway. Nevertheless, filtering as many of these as possible ensures the rest
of the pipeline runs more efficiently.

We typically include the following in our indices:

- The large and small ribosomal subunit sequences, typically from NCBI. These
  can be surprisingly difficult to track down; nevertheless, a search using the
  terms: `"ribosomal subunit"` and `organism` (including the quotes around
  `ribosomal subunit`) typically works. (click [here](https://www.ncbi.nlm.nih.gov/nuccore/?term=%22ribosomal+subunit%22+mouse) for example)

- The genomic tRNA sequences from [GtRNAdb](http://gtrnadb.ucsc.edu/). Since
  only the sequences are used, it is not very important if the annotation
  version in GtRNAdb does not match that used in the rest of the analysis. For
  example, using the tRNAs from GRCh37 and annotations from GRCh38 should not
  cause problems. Of course, if matching versions are available, it is preferred
  to use them.

- "Mt_rRNA", "Mt_tRNA" and "rRNA" genes from BioMart. In particular, we select
  those options for the "Gene type" filter. For "Attributes", we select
  "Sequences", and then specifically "Exon sequences". Additionally, including
  the "Gene type" in the header can be helpful later for identifying where reads
  in the `with-rrna` directory mapped.

Obviously, depending on the goal of the analysis, the choice of indices could
vary. For example, it may be desirable to include tRNA in the analysis or to
exclude known non-coding RNA species, such as snRNAs. These could also be
included in the "Gene type" from BioMart.

<a name="de-novo-assembly"></a>

## Using _de novo_ assembly

Often, matching RNA-seq is available for riboseq datasets. In these cases, we
highly recommend creating a _de novo_ assembly from the RNA-seq data. A variety
of algorithms are available for this task. Internally, we use a custom pipeline
which is based on the STAR aligner and StringTie assembler. Many other
options exist for transcript assembly, though. The only requirement for use with
Rp-Bp is that the assembler produces a valid GTF/GFF file (or something which can
be converted into GTF/GFF). The coordinates in the GTF/GFF file must match those in the
reference genome.

N.B. _de novo_ assembly algorithms do not typically identify coding regions (that is
what riboseq is for!); however, if they do include CDS annotations, those should
match the same start and stop codon conventions as described above.
