
# RNAxplorer
The RNAxplorer is a multitool, that offers different methods to explore RNA energy landscapes. 

## Use cases
- repellant sampling with guiding potentials on base pair level
- repellant or attractive sampling with reference structures
- retrieve local minima of secondary structures (via gradient walks in parallel)

## Description
 In default mode (or with -M RSH option) it takes an RNA sequence as input (either stdin or --sequence parameter) 
 and outputs sampled secondary RNA structures. The repellant sampling heuristic iteratively penalizes base pairs of 
 local minima of structures that have been seen too often. This results in a diverse sample set with the most
 important low free energy structures.
 
 Another important sampling method (-M RS option) is based on reference structures (--struc1
 and --struc2). This method produces structures in the vicinity of these two reference
 structures. Arbitrary many references can be added if a fasta file is used as input
 (via stdin).
 
 Often the output of sampling methods has to be coarse grained by local minima that are defined
 by a gradient walk. A parallelized gradient descent procedure can be used to retrieve
 local minima (-M RL option) of sampled structures (input via stdin).
 
 
## Example Repellant Sampling Heuristic
In order to use this method, we need a sequence or fastafile (text file) with the following content:
``` 
>sv11
GGGCACCCCCCUUCGGGGGGUCACCUCGCGUAGCUAGCUACGCGAGGGUUAAAGCGCCUUUCUCCCUCGCGUAGCUAACCACGCGAGGUGACCCCCCGAAAAGGGGGGUUUCCCA
```

In order to produce 10 sampled structures call either
```
RNAxplorer -M RSH -n 10 --sequence GGGCACCCCCCUUCGGGGGGUCACCUCGCGUAGCUAGCUACGCGAGGGUUAAAGCGCCUUUCUCCCUCGCGUAGCUAACCACGCGAGGUGACCCCCCGAAAAGGGGGGUUUCCCA
```
or
```
cat sv11.fasta | RNAxplorer -M RSH -n 10
```

The output is the the samples and the unique local minima. The columns contain the index, structure, free energy and how often this minimum was reached from a sampled structure:
``` 
GGGCACCCCCCUUCGGGGGGUCACCUCGCGUAGCUAGCUACGCGAGGGUUAAAGCGCCUUUCUCCCUCGCGUAGCUAACCACGCGAGGUGACCCCCCGAAAAGGGGGGUUUCCCA
(null)
(null)
samples so far:     10 /     10 ... done
Samples: 
(((..((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).))))))))))))))))))..))))))))...))).
(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))..))).
(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).))))))))))))))))))..))))))))).)))..
(((..((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).))))))))))))))))))))..))))))...))).
(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))..))).
(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))..))).
(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..)))))))).)))..
(((.((((((((.((((((((((((((((((...((((((((((((((..((((...))))..))))))))))))))...))))))))))))))))))...))))))))..))).
(((.(((((((((((((((((((((((((((...((((((((((((((..((((...))))..))))))))))))))...)))))))))))))))))))..))))))))..))).
(((..((((((.(((((((((((((((((((...((((((((((((((..((((...))))..))))))))))))))...)))))))))))))))))))...))))))...))).
Local minima: 
     GGGCACCCCCCUUCGGGGGGUCACCUCGCGUAGCUAGCUACGCGAGGGUUAAAGCGCCUUUCUCCCUCGCGUAGCUAACCACGCGAGGUGACCCCCCGAAAAGGGGGGUUUCCCA
   0 (((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))..))). -97.70      6
   1 (((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).))))))))))))))))))))..)))))))..))). -97.00      1
   2 (((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).))))))))))))))))))..)))))))))..))). -96.70      1
   3 (((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..)))))))).))).. -96.10      1
   4 (((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).))))))))))))))))))..))))))))).))).. -95.10      1

```

The output can also be stored in separate files:
```
cat sv11.fasta | RNAxplorer -M RSH -n 10 --lmin-file=repellant_sampling.txt
```
This creates two files with the following content:
```
repellant_sampling.txt:
     GGGCACCCCCCUUCGGGGGGUCACCUCGCGUAGCUAGCUACGCGAGGGUUAAAGCGCCUUUCUCCCUCGCGUAGCUAACCACGCGAGGUGACCCCCCGAAAAGGGGGGUUUCCCA
   0 (((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))..))). -97.70      6
   1 (((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).))))))))))))))))))))..)))))))..))). -97.00      2
   2 (((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).))))))))))))))))))..)))))))))..))). -96.70      1
   3 (((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).))))))))))))))))))..))))))))).))).. -95.10      1
```
```
repellant_sampling.samples: 
     GGGCACCCCCCUUCGGGGGGUCACCUCGCGUAGCUAGCUACGCGAGGGUUAAAGCGCCUUUCUCCCUCGCGUAGCUAACCACGCGAGGUGACCCCCCGAAAAGGGGGGUUUCCCA
(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))..))).
(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))..))).
(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))..))).
(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).))))))))))))))))))))..)))))))..))).
(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).))))))))))))))))))..)))))))))..))).
(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))..))).
(((.(((((((((((((((((((((((((((....(((((((((((((..((((...))))..)))))))))))))....)))))))))))))))))))..))))))))..))).
(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).))))))))))))))))))))..)))))))..))).
((..(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).))))))))))))))))))..)))))))))..))..
(((.(((((((((((((((((((((((((((...((((((((((((((..((((...))))..))))))))))))))...)))))))))))))))))))..))))))))..))).
```

## Example sampling with references
In order to produce structures in the vicinity of certain reference structures, you should first add the references to a fasta file:

```
>sv11
GGGCACCCCCCUUCGGGGGGUCACCUCGCGUAGCUAGCUACGCGAGGGUUAAAGCGCCUUUCUCCCUCGCGUAGCUAACCACGCGAGGUGACCCCCCGAAAAGGGGGGUUUCCCA
(((((((((((...)))))))..((((((((((....))))))))))........)))).....((((((((........)))))))).((((((((.....)))))))).....
(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))..))).

```
Then you call the following method:
```
cat sv11.fasta | RNAxplorer -M SM -e N -i 10
distortions: d_x0 = 0.3686046157, d_x1 = 0.0000000000 
1	85	-64.80	(1) ((((.((((((...))))))...((((((((((....))))))))))........)))).....((((((((........)))))))).((((((((.....)))))))).....
2	88	-65.60	(2) (((((((((((...)))))))..((((((((((....))))))))))((....)))))).....((((((((........)))))))).((((((((.....)))))))).....
2	88	-65.60	(2) (((((((((((...)))))))..((((((((((....))))))))))((....)))))).....((((((((........)))))))).((((((((.....)))))))).....
11	77	-64.00	(1) (((.(((((((...)))))))(.((((((((((....)))))))))))...(((...)))....((((((((........)))))))).((((((((.....)))))))).))).
84	2	-95.40	(1) (((.(((((((.(((((((((((((((((((...((((((((((((((..((((...))))..))))))))))))))...)))))))))))))))))))...)))))))..))).
84	4	-92.80	(1) (((.((((((((((((((((((((((((((..(.((((((((((((((..(((.....)))..)))))))))))))).)..)))))))))))))))))..)))))))))..))).
85	1	-96.70	(3) (((.(((((((((((((((((((((((((((...((((((((((((((..((((...))))..))))))))))))))...)))))))))))))))))))..))))))))..))).
85	1	-95.90	(3) (((..((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..)))))))...))).
85	1	-96.10	(3) (((.((((((((.((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).))))))))))))))))))...))))))))..))).
85	5	-94.10	(1) .((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))...)).
```
The input parameter `-M SM` selects the attractive guiding potential sampling method. With `-e N` the default distortion is used, which 
means that both reference structures have the same probability as the mfe structure. The parameter `-i 10` sets 10 iterations (1 structure per iteration will be drawn).
In order to vary the distortion in each iteration you can choose other values for `-e`. This can produce more structures on the path between two references. The
number of produced samples is then the base pair distance between both references times the number of iterations (this would be the case with `-e S`).

The output contains the computed distortions for both references. The next lines are 5 columns with the following content: base pair distance to the first reference,
base pair distance to the second reference, free energy, number of structures with the same distances to both references, the structure.
In this example you can see that samples in the base pair vicinity of both references has been constructed. In contrast to that, pure Boltzmann sampling would
produce only structures which are similar to the second reference.


## Example retrieve local minima
Local minima via gradient walks can be produced with the option `-M RL`. Input is a file with the sequence in the first line and structures in the next lines.

```
cat sv11.fasta | RNAxplorer -M RL
```

The output consists of the sequence, the index of the corresponding structure in the input file, the local minima and the energies:
```
GGGCACCCCCCUUCGGGGGGUCACCUCGCGUAGCUAGCUACGCGAGGGUUAAAGCGCCUUUCUCCCUCGCGUAGCUAACCACGCGAGGUGACCCCCCGAAAAGGGGGGUUUCCCA
2 (((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))..))). -97.70
1 (((((((((((...)))))))..((((((((((....))))))))))........)))).....((((((((........)))))))).((((((((.....))))))))..... -66.00
```




 