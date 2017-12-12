import RNA
import RNAxplorer

num_samples = 100
sequence = "GGGCACCCCCCUUCGGGGGGUCACCUCGCGUAGCUAGCUACGCGAGGGUUAAAGCGCCUUUCUCCCUCGCGUAGCUAACCACGCGAGGUGACCCCCCGAAAAGGGGGGUUUCCCA"

# init random number generator in RNAlib
RNA.init_rand()

# prepare RNAlib fold_compound
md = RNA.md()
md.uniq_ML = 1
md.compute_bpp = 0

fc = RNA.fold_compound(sequence, md, RNA.OPTION_PF)

# compute MFE structure
(ss, mfe) = fc.mfe()

# fill partition function DP matrices 
fc.pf()

for i in range(0, num_samples):
    s = fc.pbacktrack()
    print "%s" % s

# add MFE structure as first structure to repell
RNAxplorer.add_repulsion(fc, ss, -mfe / 2.)

# sample again
fc.pf()

for i in range(0, num_samples):
    s = fc.pbacktrack()
    print "%s" % s

