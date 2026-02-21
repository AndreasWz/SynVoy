#!/usr/bin/env python3
import os
import sys
import subprocess
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "bin"))
from sequence_utils import load_genome, translate, reverse_complement, write_fasta

PLOT_DIR = "results_human_ly6e_relaxed/plot_inputs_synteny_block_locus_1"
GENOMES_DIR = "results_human_ly6e_relaxed/downloaded_genomes/easy_mode_genomes"

def main():
    if not os.path.exists("alignments"):
        os.makedirs("alignments")
        
    # gene_name -> list of (species, id, prot_seq)
    extracted_seqs = defaultdict(list)
    
    for filename in os.listdir(PLOT_DIR):
        if not filename.endswith(".homology.tsv"):
            continue
            
        species = filename.split(".homology.tsv")[0]
        if species == "home_genome":
            continue
            
        # 1. Load homology map (Query_ID -> Gene_Name)
        homology_map = {}
        with open(os.path.join(PLOT_DIR, filename), "r") as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    homology_map[parts[0]] = parts[1]
                    
        # 2. Load GFF to get coordinates
        # Map Query_ID to (chrom, start, end, strand)
        gff_coords = {}
        gff_file = os.path.join(PLOT_DIR, f"{species}.gff")
        if not os.path.exists(gff_file):
            continue
            
        with open(gff_file, "r") as f:
            for line in f:
                if line.startswith("#"): continue
                parts = line.strip().split('\t')
                if len(parts) < 9: continue
                
                # Check attributes for ID
                attrs = {}
                for attr in parts[8].split(';'):
                    if '=' in attr:
                        k, v = attr.split('=', 1)
                        attrs[k] = v
                        
                q_id = attrs.get("ID")
                if q_id and q_id in homology_map:
                    # GFF is 1-based closed
                    chrom = parts[0]
                    start = int(parts[3]) - 1
                    end = int(parts[4])
                    strand = parts[6]
                    gff_coords[q_id] = (chrom, start, end, strand)
                    
        # 3. Load Genome and Extract
        
        # The species string already has .fna because the tsv is named GCF_xxx.fna.homology.tsv
        genome_file = os.path.join(GENOMES_DIR, species)
        if not os.path.exists(genome_file):
            print(f"Skipping {species}, no genome found.")
            continue
            
        print(f"Loading genome {species}...")
        genome_seqs = load_genome(genome_file)
        
        for q_id, (chrom, start, end, strand) in gff_coords.items():
            if chrom not in genome_seqs:
                continue
            dna = genome_seqs[chrom][start:end]
            if strand == '-':
                dna = reverse_complement(dna)
                
            prot = translate(dna).split('*')[0] # stop at first stop codon
            if len(prot) > 10: # Only keep reasonable length extracts
                gene_name = homology_map[q_id]
                extracted_seqs[gene_name].append((f"{species}_{q_id}", prot))
                
    # 4. Write Fasta and run MAFFT
    for gene_name, seqs in extracted_seqs.items():
        if len(seqs) < 2:
            continue
            
        fa_file = f"alignments/{gene_name}.faa"
        aln_file = f"alignments/{gene_name}_aligned.faa"
        
        fasta_entries = [(s[0], s[1]) for s in seqs]
        write_fasta(fasta_entries, fa_file)
        
        print(f"Running MAFFT for {gene_name} ({len(seqs)} sequences)...")
        with open(aln_file, "w") as out_f:
            subprocess.run(["mafft", "--auto", fa_file], stdout=out_f, stderr=subprocess.DEVNULL)
            
    print("Done! All alignments saved to alignments/")

if __name__ == "__main__":
    main()
