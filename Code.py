from Bio.SeqUtils import MeltingTemp
from Bio.Seq import Seq

# Define the target DNA sequence
target_seq = Seq("ACGACGACTACAGACGACTGACAGCAGCAGCATATACATAGAACATAGCGACAGCATAGCGACAGCATAAGCGAC")

# Design the forward primer
fwd_primer = target_seq[0:20]
fwd_tm = MeltingTemp.Tm_Wallace(fwd_primer)
print("Forward primer: %s, Tm: %.2f" % (fwd_primer, fwd_tm))

# Design the reverse primer
rev_primer = target_seq[-20:]
rev_tm = MeltingTemp.Tm_Wallace(rev_primer)
print("Reverse primer: %s, Tm: %.2f" % (rev_primer, rev_tm))

#Or
print('.........................................')


from Bio.SeqUtils import MeltingTemp
from Bio import SeqIO

# Load the target DNA sequence
target_seq = SeqIO.read('Seq.fasta', 'fasta').seq

# Define the parameters for primer design
primer_length = 20
primer_min_tm = 55
primer_max_tm = 65
primer_min_gc = 40
primer_max_gc = 60

# Design the forward primer
for i in range(len(target_seq)-primer_length):
    fwd_primer = target_seq[i:i+primer_length]
    if (MeltingTemp.Tm_Wallace(fwd_primer) >= primer_min_tm and
        MeltingTemp.Tm_Wallace(fwd_primer) <= primer_max_tm and
        100 * float(fwd_primer.count("G") + fwd_primer.count("C")) / len(fwd_primer) >= primer_min_gc and
        100 * float(fwd_primer.count("G") + fwd_primer.count("C")) / len(fwd_primer) <= primer_max_gc):
        break

# Design the reverse primer
for i in range(len(target_seq)-primer_length, 0, -1):
    rev_primer = target_seq[i:i+primer_length].reverse_complement()
    if (MeltingTemp.Tm_Wallace(rev_primer) >= primer_min_tm and
        MeltingTemp.Tm_Wallace(rev_primer) <= primer_max_tm and
        100 * float(rev_primer.count("G") + rev_primer.count("C")) / len(rev_primer) >= primer_min_gc and
        100 * float(rev_primer.count("G") + rev_primer.count("C")) / len(rev_primer) <= primer_max_gc):
        break

# Print the primers
print("Forward primer: " + fwd_primer)
print("Reverse primer: " + rev_primer)

#Tm calcluation
from Bio.SeqUtils import MeltingTemp as mt

# Define the primer sequence
Forward_primer = fwd_primer
Reverse_primer = rev_primer
# Calculate the Tm of the primer using the nearest-neighbor method
tm1 = mt.Tm_NN(Forward_primer)
tm2 = mt.Tm_NN(Reverse_primer)
# Print the Tm value
print("The Tm of the F primer is:", tm1)
print("The Tm of the R primer is:", tm2)
