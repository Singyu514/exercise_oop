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
    
class Gene(GenomicFeature):
    def __init__(self, chromosome, start, end, strand, name):
        super().__init__(chromosome, start, end, strand)
        self.name = name
        self.exons = []

    def add_exon(self, exon):
        self.exons.append(exon)

    def total_exon_length(self):
        return sum(exon.length() for exon in self.exons)

    def describe(self):
        base = super().describe()
        return f"{base.replace('Gene ', f'Gene {self.name} ')}, {len(self.exons)} exon(s)"
    
class Variant(GenomicFeature):
    def __init__(self, chromosome, start, end, strand, ref_allele, alt_allele):
        super().__init__(chromosome, start, end, strand)
        self.ref_allele = ref_allele
        self.alt_allele = alt_allele

    def variant_type(self):
        if len(self.ref_allele) == 1 and len(self.alt_allele) == 1:
            return "SNP"
        elif len(self.alt_allele) > len(self.ref_allele):
            return "insertion"
        elif len(self.alt_allele) < len(self.ref_allele):
            return "deletion"
        else:
            return "MNV"

    def describe(self):
        base = super().describe()
        return f"{base} {self.ref_allele}>{self.alt_allele} ({self.variant_type()})"

def load_data(filepath):
    genes = {}
    genes_list = []
    variants_list = []

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            record_type, chromosome, start, end, strand, field_a, field_b = fields
            start = int(start)
            end = int(end)

            if record_type == "gene":
                gene = Gene(chromosome, start, end, strand, field_a)
                genes[field_a] = gene
                genes_list.append(gene)

            elif record_type == "exon":
                exon_number = int(field_b)
                exon = Exon(chromosome, start, end, strand, exon_number)
                parent_gene = genes[field_a]
                parent_gene.add_exon(exon)

            elif record_type == "variant":
                variant = Variant(chromosome, start, end, strand, field_a, field_b)
                variants_list.append(variant)

    return genes_list, variants_list

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

    genes_list, variants_list = load_data("oop_data.tsv")

    print("\n--- Report ---")
    report_items = genes_list + variants_list
    for item in report_items:
        print(item.describe())
        if isinstance(item, Gene):
            print(f"  total exon length: {item.total_exon_length()}")

    print("\n--- Variant locations ---")
    for variant in variants_list:
        matched_genes = [g for g in genes_list if g.overlaps(variant)]

        if not matched_genes:
            print(f"{variant.describe()}: intergenic")
        else:
            for gene in matched_genes:
                in_exon = any(exon.overlaps(variant) for exon in gene.exons)
                location = "inside an exon" if in_exon else "inside gene but not in any exon"
                print(f"{variant.describe()}: inside {gene.name} ({location})")
    GenomicFeature("chr1", 5000, 1000, "+")  # should raise ValueError
