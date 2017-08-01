# codonHomologizer.pl version 1.13
Created: 6/27/2017
Last Modified: Mon Jul 31 18:05:08 2017

## WHAT IS THIS:

This script, given a set of amino acid sequences and a codon usage table, constructs the underlying genetic sequences to optimize them for crossover events.  It takes 2 or more protein sequences, aligns them to optimize codon homology (using a custom amino acid weight matrix as input to a multiple sequence alignment tool (muscle)), and then recodes the amino acid sequences to be more homologous at the DNA level.  It utilizes the supplied codon usage table to prefer codons that are more common in a particular organism.

## DETAILS:

Note, the alignment step is a multiple sequence alignment, but the sequence construction step is pair-wise.  Each sequence is produced multiple times, optimized for crossover for every other sequence.

To omit codons from incorporation in the resulting sequence, omit them from the codon usage file (or comment them out - refer to the format for the -c file below).

This script produces "Crossover Metrics" intended to be used to evaluate how good each pair of sequences is in terms of their potential to be involved in crossover events.  Here is a listing of the columns in the Crossover Metrics table:

    Source             - The file the metrics are based on.
    Pair               - Sequential/arbitrary unique pair ID.
    Aln Len            - The length of the alignment produced.
    Seqs               - Sequence IDs.
    %ID                - % identity: identicals/non-gap-length.
    %Aligned NTs       - Non-gap aln length / aln length.
    Contig Ident >= N  - Sum of contiguous identical NTs >= N.
    Rare Cdns(Usg<=N)  - Number of codons incorporated whose usage is <= N.
    Rarest Codon       - Codon:usage.  Rarest codon incorporated.
    Frame Shift bases  - Number of NTs aligned out of frame.
    Seg Aln Score[0-1] - Segment alignment score, a value between 0 and 1 (inclusive).

The Seg Aln Score is only present is the -d option is used to align for example, similar protein structural domains with one another.

To evaluate a sequence (nt alignment) produced by another tool, you can use the -e option.  All this does is add entries to the Crossover Metrics table at the end of the run.

## INPUT FORMAT: -i, --aa-file

The amino acid sequence file can be in fasta or fastq format.  There must be 2 or more sequences.  This is an unaligned sequence file.  Any alignment characters present such as gap characters ('-') or spaces will be removed.  It may be upper or lower case.  Codons for any ambiguous nucleotides will be generated as NNN.

## INPUT FORMAT: -c, --codon-file, --codon-usage-file                                                                                                                                                           

A tab- or space- delimited file of codon usage scores.  The first column is the single-character amino acid code, the second column is the codon (3 nucleotides), and the third column us the codon's usage score.  The score can be any integer or fraction.  The scores will be normalized, so the score column does not need to sum to any value.  Additional columns are optional and ignored.  Not all codons are required to be present.  All amino acids and the stop character ('*') are required to be present.  To exclude a codon from incorporation into all recoded sequences, you can either remove that row from the file or comment it out by inserting a '#' at the beginning of the line.  Empty lines are permissible.  Example:

    G	GGG	0.11	Gly
    G	GGA	0.21	Gly
    G	GGT	0.49	Gly
    G	GGC	0.19	Gly
    E	GAG	0.29	Glu
    E	GAA	0.71	Glu
    D	GAT	0.65	Asp
    D	GAC	0.35	Asp
    ...

## INPUT FORMAT: -w, --precomputed-weight-matrix                                                                           
The format of the weight matrix file is the same as is used by sequence alignment tools such as clustalw and muscle.  There is a row and column for every single-letter aminco acid code and each cell contains a score representing how good of a match the pair of amino acids is (i.e. how similar they are).  The last row and column (represented by '*') contains the minimum score.  Examples of these matrices can be found (at the time of this writing) here:

ftp://ftp.ncbi.nih.gov/blast/matrices/

Example (BLOSUM62):

       A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
    A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
    R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
    N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
    D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
    C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
    Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
    E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
    G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
    H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
    I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
    L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
    K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
    M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
    F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
    P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
    S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
    T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
    W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
    Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
    V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
    B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
    Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
    X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
    * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1
                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
## INPUT FORMAT: -e, --evaluate-nt-aln-file, --nt-file    

The aligned nucleotide sequence file can be in fasta or fastq format.  There must be 2 or more aligned sequences of the same length and they must encode a protein (i.e. the first character must be in the first 5'->3' frame).  It may be upper or lower case.
   
## INPUT FORMAT: -d, --segment-file, --domain-file                                                                                                                                 

A segment file is a tab-delimited file consisting of a series of amino acid coordinates (starting from the first amino acid as coordinate 1) in each protein.  Each row of the file is for a different sequence in the input sequence file (see -i).  The first column is the sequence ID (which must match all the IDs in the sequence file from the first non-white-space character after the defline character to the next white-space character (or the end of the line)).  E.g. a defline of ">id1 my protein" would match a value in the first column of the segment file of "id1".  See the usage for option -d for details on how the coordinates are treated.

A segment file contains coordinate pairs for each protein that are aligned with one-another independently.  There must be the same number of coordinates on each row.  An even or odd number of coordinates are acceptable, but they are treated in pairs, meaning each odd and even columned coordinate is treated as an inclusive pair such that the AAs corresponding to the pair are aligned with that segment [inclusive].  In the event that there are an odd number of coordinate columns, the coordinate in the last column is paired with the last coordinate in the protein.

Gaps between pairs are also aligned in an exclusive coordinate fashion, but if one protein contains a gap and another has no gap between pairs, e.g. 1-10 and 5-15 (where no characters in the first sequence are aligned with the first 4 characters of the second sequence), gap characters are artificially inserted.  If the first coordinate (1) or last coordinate of a protein are not included, the first and last segments are aligned such that the end coordinate is inclusive and the inner coordinate is exclusive.

Example:

    Mj	24	119	160	256	352	482	487	711
    Pf	1	100	152	269	360	544	552	769
    Sp	1	148	220	342	405	488	500	796
    Hs	38	176	230	371	419	514	516	817

Coordinates aligned:

    Mj	1-23	24-119	120-159	160-256	257-351	352-482	483-486	487-711	712-end
    Pf	*	1-100	101-151	152-269	270-359	360-544	545-551	552-769	770-end
    Sp	*	1-148	149-219	220-342	343-404	405-488	489-499	500-796	797-end
    Hs	1-37	38-176	177-229	230-371	372-418	419-514	515-515	516-817	818-end

* Sequence excluded from segment alignment and filled in with gap characters.

If a coordinate is larger than the protein provided, an error or warning will be printed.  If the second coordinate in the last pair is larger than the protein, it will be automatically adjusted to the last possible protein coordinate with a warning.  Any other coordinate out of range will generate a fatal error.

Coordinates on each row must be in sequential ascending order.

## OUTPUT FORMAT: -o, --dna-suffix, --dna-extension, --outfile, --output-file

Fasta sequence file.

## OUTPUT FORMAT: -s, --matrix-suffix, --matrix-extension                                           

The format of the weight matrix file is the same as is used by sequence alignment tools such as clustalw and muscle.  There is a row and column for every single-letter aminco acid code and each cell contains a score representing how good of a match the pair of amino acids is (i.e. how similar they are).  The last row and column (represented by '*') contains the minimum score.  Examples of these matrices can be found (at the time of this writing) here:

ftp://ftp.ncbi.nih.gov/blast/matrices/

Example (BLOSUM62):

       A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
    A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
    R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
    N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
    D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
    C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
    Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
    E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
    G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
    H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
    I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
    L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
    K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
    M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
    F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
    P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
    S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
    T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
    W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
    Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
    V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
    B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
    Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
    X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
    * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1
                                                                                                                                                                                                                                                                                                                         * OUTPUT FORMAT: --aa-aln-suffix

The format output by the tool selected by -x/-y.