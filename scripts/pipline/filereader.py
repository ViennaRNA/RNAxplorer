import re
from samplegenerator import Matrix2D

class Fasta:
    @staticmethod
    def readFile(filename):
        """
        Read a FASTA file and extract records as pairs of
        id, and sequence
        """
        res = []
        counter = 0
        is_fasta = 0
        record_id = "seq_" + `counter`
        sequence = ""
        for l in open(filename):
            l.rstrip('\n')
    
            # start actual parsing
            if l.startswith(">"):
                match = re.search('>\s*(\S+)', l)
                if match:
                    is_fasta = 1
                    record_id = match.group(1)
                    sequence = ""
            elif not l.startswith("#")  and not l.startswith(";"):
                match = re.search('([ACGUTacgutnN]+)', l)
                if match:
                    sequence += match.group(1)
    
            if not is_fasta:
                res.append((record_id, sequence))
                sequence = ""
                counter += 1
                record_id = "seq_" + `counter`
    
        return res

class xplorerFileReader:
    @staticmethod
    def readFile(filename):
        """
        Read a RNAxplorer output file and extract the references and structures.
        """
        structures = [];
        ref_struct1 = "";
        ref_struct2 = "";
        
        for l in open(filename):
            l.rstrip('\n')
    
            # start actual parsing
            match = re.search('(^[ACGUTacgutnN]+$)', l)
            if match:
                sequence = match.group(1)
            
            match = re.search('^([\(\)\.]+)\s+(\-?\d+\.\d+)$', l)
            if match:
                structure = match.group(1)
                if(ref_struct1 == ""):
                    ref_struct1 = structure
                else:
                    if (ref_struct2 == ""):
                        ref_struct2 = structure
                
            match = re.search('(^(\d+)\s+(\d+)\s+(\-?\d+\.\d+)\s+\(\d+\)\s+([\(\)\.]+)$)', l)
            if match:
                structures.append(match.group(5))
        
        result = Matrix2D(sequence, structures, ref_struct1, ref_struct2)
        return result
