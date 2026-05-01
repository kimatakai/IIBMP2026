## IIBMP2026 - Bioinformatics Programming question 1


DNA_1 = "CGTGATTCCAGCGTCTCATCAAATACGCAGCG"
DNA_2 = "CGCTGCGTATTTACGCTATTATTGGAATCACG"

RNA = "GGGGCGCAUAAAGAUGAGACGCGUUUUAGAGCUAGAAAUACGAAGUUGGGAUAAGGCUAGUGCGUUAUCAACUUGAAUUAGUGGCACCGAGUCGGUGCUUUU"

K_MER = 10


def rna2dna(seq: str) -> str:
    table = str.maketrans("AUCG", "ATCG")
    return seq.translate(table)


def reverse_dna(seq: str) -> str:
    return seq[::-1]


def complement_dna(seq: str) -> str:
    table = str.maketrans("ATCG", "TAGC")
    return seq.translate(table)



def main():
    reverse_dna_1 = reverse_dna(DNA_1)
    rna = rna2dna(RNA)
    complement_rna = complement_dna(rna2dna(RNA)) # 必要なし
    reverse_rna = reverse_dna(rna2dna(RNA))
    reverse_complement_rna = reverse_dna(complement_dna(rna2dna(RNA))) # 必要なし
    
    match_indices = []
    for rna_seq in [rna, complement_rna, reverse_rna, reverse_complement_rna]:
        for i in range(len(rna_seq) - K_MER + 1):
            kmer_fragment = rna_seq[i:i+K_MER]
            for j in range(len(reverse_dna_1) - K_MER + 1):
                if kmer_fragment == reverse_dna_1[j:j+K_MER]:
                    match_indices.append(j)
                    print(i, j, kmer_fragment, reverse_dna_1[j:j+K_MER])
    

    print(match_indices)

    if len(match_indices) == 1:
        print("Answer: ", DNA_2[match_indices[0]:match_indices[0]+K_MER])



if __name__ == "__main__":
    main()

