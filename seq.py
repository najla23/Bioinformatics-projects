#\section*{Example: Sequence Alignment of SOD1 Gene from Human and Dog}

#The following Python code demonstrates how to perform a pairwise sequence alignment between the exemplary SOD1 gene sequences of humans and dogs. The 
#sequences are aligned such that matching positions are represented with \texttt{|}.
#!/usr/bin/env python 
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio import Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

seq_human = SeqRecord(Seq("ATGGTGAGCAAGGGCGAGGAGGAGGAGAGGCGAGGAGGAGAGAGAGGAGAGAGGA"), id="Human_SOD1")
seq_dog = SeqRecord(Seq("ATGGTGAGCAAGGGCGAAGAGGAGGAGAGGGGAGGAGGAGAGAGGGGGAGAGGA"), id="Dog_SOD1")

max_len = max(len(seq_human.seq), len(seq_dog.seq))

if len(seq_human.seq) < max_len:
    seq_human.seq = seq_human.seq + '-' * (max_len - len(seq_human.seq))

if len(seq_dog.seq) < max_len:
    seq_dog.seq = seq_dog.seq + '-' * (max_len - len(seq_dog.seq))

alignment = MultipleSeqAlignment([seq_human, seq_dog])
aligner = Align.PairwiseAligner()
aligner.mode = 'global'

alignment1 = aligner.align(seq_human.seq, seq_dog.seq)
print("Alignment between Human and Dog SOD1:\n", alignment1[0])

#
#\section*{Alignment Result}
#
#The output of the alignment is as follows:
#
#\begin{verbatim}
#H-SOD1             ATGGTGAGCAAGGGCGA-GGAGGAGGAGAGGCG-AGGAGGAGAGAGAGGAG-AGAGGA
#                   |||||||||||||||||-|-|||||||||||-|-||||||||||||-||-|-||||||
#D-SOD1             ATGGTGAGCAAGGGCGAAG-AGGAGGAGAGG-GGAGGAGGAGAGAG-GG-GGAGAGGA
#\end{verbatim}


