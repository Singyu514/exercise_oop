class GenomicFeature:
    def __init__(self, chromosome, start, end, strand):
        if start <= 0 or end <= 0:
            raise ValueError("start and end must be positive integers")
        if start > end:
            raise ValueError("start must be <= end")
        if strand not in ("+", "-"):
            raise ValueError("strand must be '+' or '-'")

        self.chromosome = chromosome
        self.start = start 
        self.end = end 
        self.strand = strand 

    def length(self):
        return self.end - self.start +1
    
    def overlaps(self, other):
        if self.chromosome != other.chromosome:
            return False
        return self.start <= other.end and other.start <= self.end
    
    def describe(self):
        return f"{type(self).__name__} {self.chromosome}:{self.start}-{self.end}({self.strand})"
    
class Exon(GenomicFeature):
    def __init__(self, chromosome, start, end, strand, exon_number):
        super().__init__(chromosome, start, end, strand)
        self.exon_number = exon_number

    def describe(self):
        base = super().describe() 
        return f"{base} exon #{self.exon_number}"
    


if __name__ == "__main__":
    a = GenomicFeature("chr1", 1000, 5000, "+")
    b = GenomicFeature("chr1", 4800, 6000, "+")
    c = GenomicFeature("chr2", 1000, 5000, "+")

    print(a.describe())     # GenomicFeature chr1:1000-5000(+)
    print(a.length())        # 4001
    print(a.overlaps(b))     # True  (4800-5000 shared)
    print(a.overlaps(c))     # False (different chromosome)
   

    features = [
    GenomicFeature("chr1", 1000, 5000, "+"),
    Exon("chr1", 1000, 1200, "+", 1),
    Exon("chr1", 3000, 3300, "+", 2),
    ]
    for feature in features:
        print(feature.describe())

    GenomicFeature("chr1", 5000, 1000, "+")  # should raise ValueError

