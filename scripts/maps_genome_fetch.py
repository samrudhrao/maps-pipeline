from Bio import Entrez, SeqIO
import os

import os as _os
MAPS_HOME = "/maps" if _os.path.exists("/maps/kraken2_db") else _os.path.expanduser("~/maps")
import os, sys, subprocess, re, time

Entrez.email = "samrudh@bswhealth.org"

def fetch_genomes(organism, db_type="fungal", max_genomes=3):
    print(f"Searching NCBI for: {organism}")
    search_handle = Entrez.esearch(
        db="assembly",
        term=f"{organism}[Organism] AND latest[filter] AND complete genome[filter]",
        retmax=max_genomes)
    search_results = Entrez.read(search_handle)
    search_handle.close()
    ids = search_results["IdList"]
    if not ids:
        search_handle = Entrez.esearch(
            db="assembly",
            term=f"{organism}[Organism] AND latest[filter]",
            retmax=max_genomes)
        search_results = Entrez.read(search_handle)
        search_handle.close()
        ids = search_results["IdList"]
    if not ids:
        print(f"No genomes found for {organism}")
        return [], "0"

    fetch_handle = Entrez.efetch(db="assembly", id=",".join(ids),
                                  rettype="docsum", retmode="xml")
    fetch_results = Entrez.read(fetch_handle)
    fetch_handle.close()

    taxon_handle = Entrez.esearch(db="taxonomy",
                                   term=f"{organism}[Scientific Name]")
    taxon_results = Entrez.read(taxon_handle)
    taxon_handle.close()
    taxid = taxon_results["IdList"][0] if taxon_results["IdList"] else "0"
    print(f"  Taxid: {taxid}")

    out_dir = os.path.expanduser(f"{MAPS_HOME}/fungal_genomes")
    os.makedirs(out_dir, exist_ok=True)

    downloaded = []
    docs = fetch_results.get("DocumentSummarySet",{}).get("DocumentSummary",[])
    for doc in docs:
        try:
            acc      = str(doc.get("AssemblyAccession",""))
            ftp_path = str(doc.get("FtpPath_RefSeq",""))
            if not ftp_path or ftp_path == "na":
                ftp_path = str(doc.get("FtpPath_GenBank",""))
            if not ftp_path or ftp_path == "na": continue
            filename = f"{ftp_path.split('/')[-1]}_genomic.fna.gz"
            url      = f"{ftp_path}/{filename}"
            out_file = f"{out_dir}/{acc}_genomic.fna.gz"
            print(f"  Downloading {acc}...")
            result = subprocess.run(
                ["wget","--no-check-certificate","-q",url,"-O",out_file],
                capture_output=True)
            if result.returncode != 0 or not os.path.exists(out_file) or \
               os.path.getsize(out_file) < 1000:
                print(f"  Failed: {acc}")
                if os.path.exists(out_file): os.remove(out_file)
                continue
            subprocess.run(["gunzip","-f",out_file])
            fna_file = out_file.replace(".gz","")
            with open(fna_file) as f: content = f.read()
            content = re.sub(r'>(\S+)', f'>\\1|kraken:taxid|{taxid}', content)
            tagged = fna_file.replace(".fna","_tagged.fna")
            with open(tagged,"w") as f: f.write(content)
            os.remove(fna_file)
            downloaded.append((acc, tagged, taxid))
            print(f"  Done: {acc}")
            time.sleep(0.5)
        except Exception as e:
            print(f"  Error: {e}")
    return downloaded, taxid

def add_to_db(downloaded, db_path):
    for acc, fna, tid in downloaded:
        r = subprocess.run(
            ["kraken2-build","--add-to-library",fna,"--db",db_path],
            capture_output=True, text=True)
        status = "Added" if r.returncode==0 else "Failed"
        print(f"  {status}: {acc}")

if __name__ == "__main__":
    organism = sys.argv[1]
    db_type  = sys.argv[2] if len(sys.argv) > 2 else "fungal"
    db_path  = os.path.expanduser(f"{MAPS_HOME}/kraken2_db_fungal_refseq")
    result, taxid = fetch_genomes(organism, db_type)
    if result:
        print(f"\nAdding {len(result)} genomes to database...")
        add_to_db(result, db_path)
        print(f"Done — {len(result)} genomes added for {organism}")
    else:
        print(f"No genomes downloaded for {organism}")
