# AliHassan Hashim
# 8/4/25
# Professor Larkin
#
# Due: 8/6/25
#
# Final - DNA/RNA Strand Analysis
#
# This program reads in lines from DnaRna_Library.txt and determines mutations, 
# DNA/RNA sequences, and noise regions.
# The results are formatted and printed.

class Base:
    def __init__(self, base):
        self.base = base

    # This function determines if a base is valid depending on the type of strand it is (RNA/DNA/invalid)
    def is_valid_for(self, type):
        options = []
        if type == "DNA":
            options = ['A', 'T', 'C', 'G']
        elif type == "RNA":
            options = ['A', 'U', 'C', 'G']
        else:
            raise ValueError(f"{type} is an invalid type of strand")

        if self.base not in options:
            return False
        return True

    # This function checks if two bases are equal for comparison
    def __eq__(self, other_base):
        return self.base == other_base.base

    # Returns info as printable string
    def __str__(self):
        return self.base

    def get_base(self):
        return self.base


class Codon:
    def __init__(self, bases):
        self.bases = bases      # List of base objects

    # This function ensures the codon is valid (length and bases)
    def is_valid(self):  
        if len(self.bases) != 3:
            raise ValueError("Error, the length of each codon must be 3.")
        
        for base in self.bases:
            if base not in ['A', 'T', 'C', 'G', 'U']:
                raise ValueError(f"Error: The base {base} is invalid")
            
        return True
    
    # This function creates a string to be printed for the codon
    def to_str(self):
        codon_str = ""
        for base in self.bases:
            codon_str += str(base)
        return codon_str
    
    # This function allows for comparison of two codons
    def __eq__(self, other):
        if len(self.bases) != len(other.bases):
            return False
        for index in range(0, len(self.bases)):
            if self.bases[index] != other.bases[index]:
                return False
        return True
    
    # Returns info as printable string
    def __str__(self):
        return self.to_str()
    

class ProteinRegion:
    def __init__(self, codons):
        self.codons = codons

    # Converts class info to printable string
    def to_str(self):
        segment_str = ""
        for codon in self.codons:
            segment_str += codon.to_str()
        return segment_str


class Segment:
    def __init__(self, strand):
        self.strand = strand
        self.bases = []

        for letter in self.strand:
            self.bases.append(Base(letter))

    # These functions will be overridden in child classes
    def is_valid(self):
        pass

    # Returns info as printable string
    def to_str(self):
        pass

    # This function formats the segments for printing
    def format_codons(self):
        codons = self.get_codons()
        out_str = ""
        for codon in codons:
            out_str += str(codon) + " "
        return out_str.strip()
    
    # This function returns the strand as codons
    def get_codons(self):
        codons = []
        for index in range(0, len(self.bases), 3):      # Separate every 3 bases
            codons.append(Codon(self.bases[index:index+3]))

        return codons

    # This function extracts the index values of RNA and DNA protein regions and returns them as a list
    def find_protein_regions(self):
        regions_index = []
        codons = self.get_codons()

        for index in range(0, len(codons)):
            codon = codons[index]
            start_index = 0

            if self.start_codon(codon) != "":
                type = self.start_codon(codon)  # DNA or RNA strand, allows proper stop codon to be searched for

                start_index = index

                index_2 = index
                codon_2 = codon
                
                while index_2 < len(codons) and not self.stop_codon(codon_2, type): # Search for end of region
                    index_2 += 1
                    codon_2 = codons[index_2]
                
                stop_index = index_2
                regions_index.append([start_index, stop_index])
        return regions_index

    # This function finds the start codon in a string and determines if it is DNA or RNA
    def start_codon(self, codon):
        if codon.to_str() == "ATG":
            return "DNA"
        if codon.to_str() == "AUG":
            return "RNA"
        return ""
    
    # This function finds the stop codon in a string depending on if it is DNA or RNA
    def stop_codon(self, codon, type):
        if type == "DNA":
            if codon.to_str() in ["TAA", "TAG", "TGA"]:
                return True
        if type == "RNA":
            if codon.to_str() in ["UAA", "UAG", "UGA"]:
                return True
        return False
    
    # Returns info as printable string
    def __str__(self):
        return self.to_str()


class DNA(Segment):
    # This function ensures the segment is a valid DNA region
    def is_valid(self):
        codons = self.get_codons()
        if self.start_codon(codons[0]) != "DNA":
            return False
        elif not self.stop_codon(codons[-1], "DNA"):
            return False
        
        for base in self.bases:
            if not base.is_valid_for("DNA"):
                return False

        return True
    
    # Returns info as printable string
    def to_str(self):
        return self.format_codons()

class RNA(Segment):
    # This function ensures the segment is a valid RNA region
    def is_valid(self):
        codons = self.get_codons()
        if self.start_codon(codons[0]) != "RNA":
            return False
        elif not self.stop_codon(codons[-1], "RNA"):
            return False
        
        for base in self.bases:
            if not base.is_valid_for("RNA"):
                return False
            
        return True
    
    # Returns info as printable string
    def to_str(self):
        return self.format_codons()

class Noise(Segment):
    # Returns info as printable string
    def to_str(self):
        return "[NOISE]" + self.format_codons() + "[/NOISE]"

class Mutation(Segment):
    # This function ensures the segment is a valid mutated region
    def is_valid(self):
        if "T" in self.strand and "U" in self.strand:
            return True
        return False

    # Returns info as printable string
    def to_str(self):
        return "[MUT]" + self.format_codons() + "[/MUT]"
    

class StrandFactory:
    def __init__(self, str):
        self.str = str

    # Extracts/slices the segments of DNA/RNA protein regions from the strand, isolates noise
    def get_segments(self):
        segments = self.slice()
        classified_str = []
        for segment in segments:
            str_class = self.classify(segment)
            classified_str.append(str_class)

        return classified_str

    # This function classifies the type of segment: DNA, RNA, Noise, or Mutation
    def classify(self, segment):
        dna_obj = DNA(segment)
        rna_obj = RNA(segment)
        mut_obj = Mutation(segment)
        noise_obj = Noise(segment)

        if dna_obj.is_valid():
            return dna_obj
        elif rna_obj.is_valid():
            return rna_obj
        elif mut_obj.is_valid():
            return mut_obj
        else:
            return noise_obj

    # This function separates the segments from the main strand
    def slice(self):
        main_segment = Segment(self.str)
        protein_regions = main_segment.find_protein_regions()
        protein_strands = []
        for values in protein_regions:
            protein_strands.append(self.str[values[0] * 3 : (values[1] + 1) * 3])   # Separate protein regions
        
        new_str = self.str
        index = 0

        for strand in protein_strands:
            start_ind = new_str.find(strand, index)
            stop_ind = start_ind + len(strand)
            
            new_str = new_str[0:start_ind] + " " + new_str[start_ind:stop_ind] + " " + new_str[stop_ind:]
            # Separates segments using whitespace

            index += 1
        new_str = new_str.strip()

        split_str = new_str.split()
        
        return split_str
    
class Strand:
    def __init__(self, strand):
        self.strand = strand

        str_fact = StrandFactory(self.strand)
        self.sliced_segments = str_fact.get_segments()

    def __str__(self):
        output = ""
        for segment in self.sliced_segments:
            output += segment.to_str() + " "
        return output
    
    def get_segments(self):
        return self.sliced_segments
    
    def render(self):
        return str(self).replace(" ", "")
    
    def to_str(self):
        return str(self)


class StrandAnalyzer:
    def __init__(self, strand_obj):
        self.strand_obj = strand_obj
        self.segments = strand_obj.get_segments()

    # This function extracts codons from all valid segments (DNA or RNA only)
    def get_valid_codons(self):
        valid_codons = []

        for segment in self.segments:
            if isinstance(segment, DNA) or isinstance(segment, RNA):
                seg_codons = segment.get_codons()

                for codon in seg_codons:
                    valid_codons.append(codon)

        return valid_codons
                

    # This function extracts invalid codons from valid DNA/RNA sequences
    def get_mutations(self):
        mutated_codons = []

        for segment in self.segments:
            if isinstance(segment, Mutation):
                seg_codons = segment.get_codons()

                for codon in seg_codons:
                    mutated_codons.append(codon)

        return mutated_codons

    # This function extracts codons from non DNA/RNA sequences
    def get_noise(self):
        noise_codons = []

        for segment in self.segments:
            if isinstance(segment, Noise):
                seg_codons = segment.get_codons()
                
                for codon in seg_codons:
                    noise_codons.append(codon)

        return noise_codons

# This function prints the different types of segment codons found in the main strand
def print_codons(strand_an_call):
    valid_codons = strand_an_call.get_valid_codons()
    mutated_codons = strand_an_call.get_mutations()
    noise_codons = strand_an_call.get_noise()

    valid_str = ""

    for valid_codon in valid_codons:
        valid_str += valid_codon.to_str() + " "

    mutated_str = ""

    for mutated_codon in mutated_codons:
        mutated_str += mutated_codon.to_str() + " "

    noise_str = ""

    for noise_codon in noise_codons:
        noise_str += noise_codon.to_str() + " "

    print(f"\nValid Codons:\n{valid_str}")
    print(f"\nMutated Codons:\n{mutated_str}")
    print(f"\nNoise Codons:\n{noise_str}")


# This function takes DNA info from the DNA_Library.txt and combines it all into one string
def extract_info(filename):
    strand = ""

    try:
        with open(filename, 'r') as file:
            for line in file:
                strand += line.strip()  # Combines into one big strand of DNA
        file.close()
        strand = strand.strip()

    except FileNotFoundError:       # Handles file not being found
        raise ValueError(f"Error: File {filename} not found")

    return strand

# The main method to run the entire program
def main():
    strand_in = extract_info("DnaRna_Library.txt")
    strand_call = Strand(strand_in)
    print("Complete Strand:\n" + strand_call.to_str())
    strand_an_call = StrandAnalyzer(strand_call)
    print_codons(strand_an_call)

if __name__ == "__main__":
    main()