
import os
import sys
import ftplib
import urllib.request
import re

GENOMES = [
    ("Homo_sapiens", "GCF_000001405.40"),
    ("Anolis_carolinensis", "GCF_000090745.1"),
    ("Naja_naja", "GCA_009733165.1"),
    ("Pseudonaja_textilis", "GCA_900518735.1"),
    ("Protobothrops_mucrosquamatus", "GCA_001527695.3"),
    ("Python_bivittatus", "GCA_000186305.2"),
    ("Bungarus_multicinctus", "GCA_023653725.1")
]

OUT_DIR = "genomes"

def get_ftp_path(accession):
    """Constructs the FTP path from accession."""
    prefix, numbers = accession.split('_')
    version_split = numbers.split('.')
    base_nums = version_split[0]
    
    # GCF/000/001/405/
    part1 = base_nums[0:3]
    part2 = base_nums[3:6]
    part3 = base_nums[6:9]
    
    return f"/genomes/all/{prefix}/{part1}/{part2}/{part3}"

def download_genome(name, accession):
    print(f"[{name}] Resolving {accession}...")
    ftp_path = get_ftp_path(accession)
    
    try:
        ftp = ftplib.FTP('ftp.ncbi.nlm.nih.gov')
        ftp.login()
        ftp.cwd(ftp_path)
        
        # List files to find the specific folder
        try:
            files = ftp.nlst()
        except ftplib.error_perm as resp:
            if str(resp) == "550 No files found":
                print("No files in this directory")
                files = []
            else:
                raise
        
        # Find directory starting with accession
        target_dir = None
        for fname in files:
            # NLST returns filenames/dirnames in current dir
            if fname.startswith(accession):
                target_dir = fname
                break
        
        if not target_dir:
            print(f"Error: Could not find directory for {accession} in {ftp_path}")
            return
            
        ftp.cwd(target_dir)
        
        # Find genomic.fna.gz
        try:
            files = ftp.nlst()
        except ftplib.error_perm:
            print(f"Error listing files in {target_dir}")
            return
        
        target_file = None
        for fname in files:
            # NLST returns relative names
            if fname.endswith("_genomic.fna.gz") and not "from_genomic" in fname:
                target_file = fname
                break
                
        if not target_file:
            print(f"Error: Could not find genomic fasta in {target_dir}")
            return
            
        # Download
        download_url = f"https://ftp.ncbi.nlm.nih.gov{ftp_path}/{target_dir}/{target_file}"
        output_path = os.path.join(OUT_DIR, f"{name}.fna.gz")
        
        if os.path.exists(output_path):
            print(f"[{name}] Already exists.")
            return

        print(f"[{name}] Downloading {download_url}...")
        urllib.request.urlretrieve(download_url, output_path)
        print(f"[{name}] Done.")
        
        # Gunzip? Pipeline handles gz?
        # main.nf handles .gz usually if configured, but let's gunzip to be safe or check config.
        # Actually nextflow usually handles gz seamlessly.
        
    except Exception as e:
        print(f"Error processing {name}: {e}")

if __name__ == "__main__":
    if not os.path.exists(OUT_DIR):
        os.makedirs(OUT_DIR)
        
    for name, acc in GENOMES:
        download_genome(name, acc)
