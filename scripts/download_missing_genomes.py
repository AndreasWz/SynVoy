#!/usr/bin/env python3
import os
import ftplib
import urllib.request
import subprocess

GENOMES = [
    ("Naja_naja", "GCA_009733165.1"),
    ("Python_bivittatus", "GCA_000186305.2"),
    ("Bungarus_multicinctus", "GCA_023653725.1"),
    ("Thamnophis_sirtalis", "GCA_001077635.1"),
    ("Myotis_lucifugus", "GCA_000147115.1") # Optional, just fetching snakes
]

OUT_DIR = "combined_test_genomes"

def get_ftp_path(accession):
    prefix, numbers = accession.split('_')
    version_split = numbers.split('.')
    base_nums = version_split[0]
    return f"/genomes/all/{prefix}/{base_nums[0:3]}/{base_nums[3:6]}/{base_nums[6:9]}"

def download_genome(name, accession):
    print(f"Resolving {accession}...")
    ftp_path = get_ftp_path(accession)
    
    try:
        ftp = ftplib.FTP('ftp.ncbi.nlm.nih.gov')
        ftp.login()
        ftp.cwd(ftp_path)
        
        target_dir = [d for d in ftp.nlst() if d.startswith(accession)]
        if not target_dir:
            return
        target_dir = target_dir[0]
        ftp.cwd(target_dir)
        
        files = ftp.nlst()
        fna_file = next((f for f in files if f.endswith("_genomic.fna.gz") and "from_genomic" not in f), None)
        gff_file = next((f for f in files if f.endswith("_genomic.gff.gz")), None)
        
        for file_to_dl, ext in [(fna_file, "fna"), (gff_file, "gff")]:
            if not file_to_dl: continue
            
            download_url = f"https://ftp.ncbi.nlm.nih.gov{ftp_path}/{target_dir}/{file_to_dl}"
            base_out = os.path.join(OUT_DIR, f"{accession}.{ext}")  # naming them by accession to match others
            gz_out = f"{base_out}.gz"
            
            if not os.path.exists(base_out):
                print(f"Downloading {download_url} to {gz_out}...")
                urllib.request.urlretrieve(download_url, gz_out)
                subprocess.run(["gunzip", gz_out])
                print(f"Extracted to {base_out}")
            else:
                print(f"{base_out} already exists.")
                
    except Exception as e:
        print(f"Failed {accession}: {e}")

if __name__ == "__main__":
    for name, acc in GENOMES:
        download_genome(name, acc)
