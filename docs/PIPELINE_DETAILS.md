# SynTerra Pipeline - Detaillierte Technische Dokumentation

## Überblick

SynTerra ist eine Nextflow-basierte Pipeline zur **Syntenie-Analyse von Genen über mehrere Genome hinweg**. Die Pipeline lokalisiert ein "Gene of Interest" (GOI) in einem Referenzgenom, identifiziert flankierende Gene, und sucht dann nach orthologen Regionen in verwandten Genomen.

---

## Kernkonzept: Was ist Syntenie?

Syntenie bedeutet, dass Gene in der gleichen Reihenfolge auf dem Chromosom liegen. Wenn ein GOI in Genom A von den Genen X-Y-Z flankiert wird, und in Genom B findet man ebenfalls X-Y-Z in dieser Reihenfolge, spricht das für:
- Gemeinsame evolutionäre Herkunft
- Konservierte genomische Organisation
- Höhere Wahrscheinlichkeit echter Orthologie

---

## Pipeline-Ablauf im Detail

### Phase 1: Eingabe und Vorbereitung

#### 1.1 FETCH_QUERY
**Zweck**: Lädt die Query-Sequenz (das GOI)

```
Input:  params.gene (Pfad zu FASTA oder Accession-Nummer)
Output: gene_query.fasta
```

Wenn eine Accession-Nummer angegeben wird (z.B. `NP_001011614.1`), wird die Sequenz von NCBI heruntergeladen. Sonst wird die lokale FASTA-Datei verwendet.

#### 1.2 FETCH_RELATED_GENOMES
**Zweck**: Lädt verwandte Genome von NCBI basierend auf Taxonomie

```python
# bin/fetch_related_genomes.py
# Sucht in NCBI Assembly nach Genomen der angegebenen taxonomischen Gruppe
# Filter: nur RefSeq, nur "Complete Genome" oder "Chromosome"
```

**Parameter**:
- `params.taxon`: Taxonomische Gruppe (z.B. "Formicidae" für Ameisen)
- `params.max_genomes`: Maximale Anzahl zu ladender Genome

#### 1.3 STAGE_GENOMES
**Zweck**: Kombiniert lokale und heruntergeladene Genome

```
Input:  - Lokale Genome aus params.genomes_dir
        - Heruntergeladene Genome von FETCH_RELATED_GENOMES
Output: Tuple (genome_id, genome.fna, genome.gff oder null)
```

---

### Phase 2: GOI-Lokalisierung im Home-Genom

#### 2.1 LOCATE_GENE
**Zweck**: Findet die Position des GOI im Home-Genom

```bash
# Zwei redundante Suchstrategien für maximale Sensitivität:

# 1. MMseqs2 mit search-type 2 (= tblastn-Modus)
mmseqs easy-search query.faa home_genome.fna hits_mmseqs.m8 tmp \
    --search-type 2 \
    --min-seq-id 0.3 \
    -e 1e-5

# 2. BLAST tblastn als Backup
tblastn -query query.faa -subject home_genome.fna -outfmt 6
```

**Wichtig**: Der `--search-type 2` Parameter bei MMseqs2 bedeutet:
- Query = Protein
- Target = Nukleotid (6-Frame-Translation)
- Entspricht BLAST's tblastn

**Output**: `home_gene_location.bed` mit den Koordinaten des GOI

```
# Format: Chromosom  Start  Ende  Name  Score  Strand
NC_037638.1    1234567    1235890    melittin    95.5    +
```

---

### Phase 3: Home-Genom Proteom-Vorbereitung

#### 3.1 PREPARE_HOME_PROTEOME
**Zweck**: Erstellt ein Proteom des Home-Genoms

**Fall A: GFF-Datei vorhanden**
```bash
# Extrahiere Protein-Sequenzen aus dem Genom basierend auf GFF-Annotation
python gff_to_faa.py --genome home.fna --gff home.gff --output home_proteins.faa
```

**Fall B: Keine GFF-Datei (nur Genom)**
```bash
# Verwende Prodigal für ab-initio Genvorhersage
prodigal -i home.fna -a home_proteins.faa -f gff -o home_predicted.gff -p meta
```

**Output**:
- `home_proteins.faa`: Alle vorhergesagten/annotierten Proteine
- `home_predicted.gff`: GFF-Datei (falls durch Prodigal generiert)

---

### Phase 4: Flankierende Gene extrahieren

#### 4.1 EXTRACT_FLANKING_GENES
**Zweck**: Identifiziert Gene um das GOI herum

```python
# bin/extract_flanking_genes.py

# 1. Lade GOI-Position aus BED-Datei
goi_chrom = "NC_037638.1"
goi_start = 1234567
goi_end = 1235890

# 2. Parse GFF-Datei für alle Gene auf demselben Chromosom
genes_on_chrom = [g for g in all_genes if g.chrom == goi_chrom]

# 3. Sortiere nach Position
genes_sorted = sorted(genes_on_chrom, key=lambda g: g.start)

# 4. Finde GOI-Index und extrahiere Nachbarn
goi_idx = find_closest_gene(genes_sorted, goi_start, goi_end)
flanking_genes = genes_sorted[goi_idx - N : goi_idx + N + 1]  # N = params.flanking_genes

# 5. Extrahiere Protein-Sequenzen der flankierenden Gene
```

**Parameter**: `params.flanking_genes` (default: 5) - Anzahl Gene auf jeder Seite

**Output**: 
- `flanking_proteins.faa`: Proteinsequenzen der flankierenden Gene
- `flanking_genes.tsv`: Metadaten (Position, Strand, Produkt-Name)

---

### Phase 5: Iterative Suche in Zielgenomen

#### 5.1 ITERATIVE_SEARCH (Kernprozess)
**Zweck**: Sucht orthologe Regionen in jedem Zielgenom

```python
# bin/iterative_search_runner.py - Hauptlogik

def search_genome(target_genome, query_proteins, home_proteins):
    """
    Für jedes Zielgenom:
    1. Suche mit Query-Proteinen
    2. Validiere durch Reciprocal Best Hit
    3. Annotiere Treffer mit Miniprot
    """
```

##### Schritt 5.1.1: Initiale Suche mit MMseqs2

```bash
# Proteine gegen Genom suchen (tblastn-Modus)
mmseqs easy-search flanking_proteins.faa target_genome.fna hits.m8 tmp \
    --search-type 2 \
    -e 1e-5 \
    --min-seq-id 0.3
```

##### Schritt 5.1.2: Reciprocal Best Hit (RBH) Validierung

```python
def validate_rbh(hit, target_genome, home_proteome):
    """
    Validiert ob ein Treffer ein echter Ortholog ist.
    
    1. Extrahiere die gefundene Sequenz aus dem Zielgenom
    2. Suche diese Sequenz gegen das Home-Proteom
    3. Prüfe: Ist der beste Treffer das ursprüngliche Query-Protein?
    
    Wenn ja → Reciprocal Best Hit → wahrscheinlich echter Ortholog
    Wenn nein → Möglicherweise Paralog oder false positive
    """
    
    # Vorwärtssuche: Query → Target (bereits erfolgt)
    forward_hit = hit  # z.B. "melittin" findet "region_12345" in Target
    
    # Rückwärtssuche: Target-Region → Home-Proteom
    target_sequence = extract_sequence(target_genome, hit.start, hit.end)
    reverse_hits = mmseqs_search(target_sequence, home_proteome)
    
    # Prüfung
    best_reverse_hit = reverse_hits[0]
    if best_reverse_hit.target == hit.query:
        return True  # RBH bestätigt!
    else:
        return False  # Kein RBH
```

##### Schritt 5.1.3: Miniprot-Annotation

**Was ist Miniprot?**
Miniprot ist ein Protein-zu-Genom-Aligner, der:
- Intron-Exon-Strukturen erkennt
- Splice-Sites findet
- Proteine auf genomische DNA mappt

```python
def run_miniprot(target_fasta, query_protein, output_paf):
    """
    Führt Miniprot aus und parst die Ergebnisse.
    
    Miniprot-Befehl:
    miniprot -I --gff target.fna query.faa
    
    -I: Verhindert Alignments über Contigs hinweg
    --gff: Gibt GFF-Format aus (enthält CDS-Koordinaten)
    """
    
    cmd = ["miniprot", "-I", "--gff", target_fasta, query_protein]
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    # Parse GFF-Output
    hits = []
    current_hit = None
    
    for line in result.stdout.split('\n'):
        if line.startswith('#'):
            continue
        
        fields = line.split('\t')
        feature_type = fields[2]  # mRNA, CDS, etc.
        
        if feature_type == 'mRNA':
            # Neue mRNA = neuer Hit
            current_hit = {
                'chrom': fields[0],
                'start': int(fields[3]),
                'end': int(fields[4]),
                'strand': fields[6],
                'score': parse_score(fields[8]),
                'cds_parts': []  # Hier kommen die Exons rein
            }
            hits.append(current_hit)
            
        elif feature_type == 'CDS':
            # CDS = ein Exon des aktuellen Hits
            cds_start = int(fields[3])
            cds_end = int(fields[4])
            current_hit['cds_parts'].append((cds_start, cds_end))
    
    return hits
```

##### Schritt 5.1.4: CDS-Extraktion (Exon-Level)

```python
def extract_cds_sequence(genome_seq, hit):
    """
    Extrahiert die kodierender Sequenz auf Exon-Level.
    
    Wichtig: Ein Gen kann aus mehreren Exons bestehen!
    
    Beispiel:
    Gen auf + Strang:
    Exon1: 1000-1200 (200 bp)
    Exon2: 1500-1800 (300 bp)  
    Exon3: 2000-2150 (150 bp)
    
    CDS = Exon1 + Exon2 + Exon3 = 650 bp
    """
    
    # Sortiere Exons nach Position
    sorted_cds = sorted(hit['cds_parts'], key=lambda x: x[0])
    
    # Extrahiere und konkateniere DNA-Sequenzen
    cds_sequence = ""
    for start, end in sorted_cds:
        exon_seq = genome_seq[start-1:end]  # 1-based zu 0-based
        cds_sequence += exon_seq
    
    # Strand-Handling
    if hit['strand'] == '-':
        # Minus-Strang: Reverse Complement
        cds_sequence = reverse_complement(cds_sequence)
    
    # Translation zu Protein
    protein = translate(cds_sequence)
    
    return cds_sequence, protein
```

**Warum ist das wichtig?**

Ohne Exon-Level Extraktion würde man:
```
Falsch: Genom[start:end] → enthält Introns!
        ATGCCC...INTRON...GGGAAA...INTRON...TTTAAA

Richtig: Exon1 + Exon2 + Exon3 → nur kodierende Sequenz
         ATGCCC + GGGAAA + TTTAAA
```

##### Schritt 5.1.5: Clustering und Region-Definition

```python
def cluster_hits(hits, max_gap=50000):
    """
    Gruppiert nahe beieinander liegende Hits zu genomischen Regionen.
    
    Wenn mehrere flankierende Gene in der Nähe gefunden werden,
    gehören sie wahrscheinlich zur selben syntenischen Region.
    """
    
    # Sortiere nach Chromosom und Position
    sorted_hits = sorted(hits, key=lambda h: (h['chrom'], h['start']))
    
    clusters = []
    current_cluster = [sorted_hits[0]]
    
    for hit in sorted_hits[1:]:
        last_hit = current_cluster[-1]
        
        # Gleicher Chromosom und nahe genug?
        if (hit['chrom'] == last_hit['chrom'] and 
            hit['start'] - last_hit['end'] < max_gap):
            current_cluster.append(hit)
        else:
            # Neuer Cluster
            clusters.append(current_cluster)
            current_cluster = [hit]
    
    clusters.append(current_cluster)
    return clusters
```

**Output von ITERATIVE_SEARCH**:
- `{genome_id}_hits.bed`: Alle gefundenen Hits mit Koordinaten
- `{genome_id}_proteins.faa`: Extrahierte Proteinsequenzen
- `{genome_id}_cds.fna`: Extrahierte CDS-Nukleotidsequenzen
- `{genome_id}_annotations.gff`: GFF-Annotation der Treffer

---

### Phase 6: Augmented Search (Optional)

#### 6.1 AUGMENTED_SEARCH
**Zweck**: Erweiterte Suche mit neu gefundenen Homologen

```python
# Idee: Verwende gefundene Orthologe als zusätzliche Queries

# Runde 1: Suche mit Original-Query
# → Findet: Ortholog_A, Ortholog_B

# Runde 2: Suche mit Original + Ortholog_A + Ortholog_B
# → Findet möglicherweise: Ortholog_C (der dem Original zu unähnlich war)

# Das erhöht die Sensitivität bei divergenten Sequenzen
```

---

### Phase 7: Analyse und Phylogenie

#### 7.1 CLUSTER_REGIONS
**Zweck**: Gruppiert alle gefundenen Regionen

```python
# bin/cluster_grs.py
# Fasst redundante Treffer zusammen
# Identifiziert Paraloge vs. Orthologe
```

#### 7.2 COMPUTE_TREE
**Zweck**: Berechnet phylogenetische Bäume

```bash
# 1. Multiple Sequence Alignment mit MAFFT
mafft --auto all_proteins.faa > aligned.faa

# 2. Phylogenetischer Baum mit FastTree
fasttree aligned.faa > tree.nwk
```

#### 7.3 PHYLO_SORT
**Zweck**: Sortiert Ergebnisse nach Phylogenie

```python
# bin/phylo_sort.py
# Ordnet Treffer basierend auf dem phylogenetischen Baum
# Identifiziert Clades von Orthologen
```

---

### Phase 8: Visualisierung

#### 8.1 PLOT_SYNTENY
**Zweck**: Erstellt interaktive Syntenie-Plots

```python
# bin/plot_synteny.py

# Input:
# - Home-GFF: Zeigt flankierende Gene im Referenzgenom
# - Target-BED-Dateien: Zeigen gefundene Regionen in Zielgenomen
# - Phylogenetischer Baum: Für Sortierung

# Output:
# - synteny_plot.html: Interaktiver Plotly-Plot
# - synteny_plot.png: Statisches Bild
```

**Der Plot zeigt**:
- Horizontale Tracks für jedes Genom
- Pfeile für Gene (Richtung = Strand)
- Farben für orthologe Gruppen
- Verbindungslinien zwischen Orthologen

#### 8.2 ANNOTATE_STRUCTURE
**Zweck**: Annotiert Proteindomänen

```bash
# Suche gegen Pfam/InterPro Datenbanken
hmmscan --domtblout domains.txt Pfam-A.hmm proteins.faa
```

---

### Phase 9: Report-Generierung

#### 9.1 GENERATE_REPORT
**Zweck**: Erstellt den finalen HTML-Report

```python
# bin/generate_report.py

# Kombiniert alle Ergebnisse:
# - Statistiken (Anzahl Genome, gefundene Gene)
# - Phylogenetischer Baum (interaktiv)
# - Syntenie-Plot
# - Sequenz-Alignments
# - Domänen-Annotationen
```

---

## Datenfluss-Diagramm

```
                    ┌─────────────────┐
                    │   INPUT         │
                    │  - GOI (FASTA)  │
                    │  - Home Genom   │
                    │  - Taxon        │
                    └────────┬────────┘
                             │
              ┌──────────────┴──────────────┐
              ▼                             ▼
    ┌─────────────────┐           ┌─────────────────┐
    │  FETCH_QUERY    │           │ FETCH_RELATED   │
    │  (Lade GOI)     │           │    GENOMES      │
    └────────┬────────┘           └────────┬────────┘
             │                             │
             ▼                             ▼
    ┌─────────────────┐           ┌─────────────────┐
    │  LOCATE_GENE    │           │ STAGE_GENOMES   │
    │ (Finde GOI im   │           │ (Sammle alle    │
    │  Home-Genom)    │           │  Genome)        │
    └────────┬────────┘           └────────┬────────┘
             │                             │
             ▼                             │
    ┌─────────────────┐                    │
    │ PREPARE_HOME    │                    │
    │   PROTEOME      │                    │
    └────────┬────────┘                    │
             │                             │
             ▼                             │
    ┌─────────────────┐                    │
    │EXTRACT_FLANKING │                    │
    │     GENES       │                    │
    └────────┬────────┘                    │
             │                             │
             └──────────────┬──────────────┘
                            ▼
                  ┌─────────────────┐
                  │ITERATIVE_SEARCH │
                  │                 │
                  │ Für jedes Genom:│
                  │ 1. MMseqs2      │
                  │ 2. RBH-Check    │
                  │ 3. Miniprot     │
                  │ 4. CDS-Extract  │
                  └────────┬────────┘
                           │
              ┌────────────┴────────────┐
              ▼                         ▼
    ┌─────────────────┐       ┌─────────────────┐
    │ CLUSTER_REGIONS │       │ AUGMENTED_SEARCH│
    └────────┬────────┘       └────────┬────────┘
             │                         │
             └────────────┬────────────┘
                          ▼
                ┌─────────────────┐
                │  COMPUTE_TREE   │
                │  (MAFFT +       │
                │   FastTree)     │
                └────────┬────────┘
                         │
                         ▼
                ┌─────────────────┐
                │   PHYLO_SORT    │
                └────────┬────────┘
                         │
              ┌──────────┴──────────┐
              ▼                     ▼
    ┌─────────────────┐   ┌─────────────────┐
    │  PLOT_SYNTENY   │   │ANNOTATE_STRUCT  │
    └────────┬────────┘   └────────┬────────┘
             │                     │
             └──────────┬──────────┘
                        ▼
              ┌─────────────────┐
              │ GENERATE_REPORT │
              │   (HTML)        │
              └─────────────────┘
```

---

## Wichtige Algorithmen im Detail

### Reciprocal Best Hit (RBH)

```
Genom A (Home)              Genom B (Target)
     │                            │
     │  Protein X ───────────────►│ Region Y (bester Hit)
     │                            │
     │◄─────────────────── Region Y sucht zurück
     │                            │
     │  Bester Treffer = Protein X?
     │        │
     │        ├── JA → RBH bestätigt → Orthologe!
     │        └── NEIN → Möglicherweise Paralog
```

### Miniprot Splice-Alignment

```
Query Protein:    MKFLILLFNILVSAP...

Genom:            ═══════════╗         ╔════════════╗         ╔═══════════
                  ║  Exon 1  ║─intron──║   Exon 2   ║─intron──║  Exon 3
                  ╚══════════╝         ╚════════════╝         ╚═══════════

Miniprot Output:
  mRNA  chr1  1000  3000  +  (Gesamtregion)
  CDS   chr1  1000  1200  +  (Exon 1)
  CDS   chr1  1500  1800  +  (Exon 2)
  CDS   chr1  2700  3000  +  (Exon 3)
```

---

## Konfigurationsparameter

| Parameter | Default | Beschreibung |
|-----------|---------|--------------|
| `gene` | required | GOI als FASTA oder Accession |
| `home_genome` | required | Referenzgenom (FASTA) |
| `home_gff` | null | GFF-Annotation (optional) |
| `taxon` | required | Taxonomische Gruppe für Suche |
| `max_genomes` | 10 | Maximale Anzahl Genome |
| `flanking_genes` | 5 | Gene auf jeder Seite des GOI |
| `evalue` | 1e-5 | E-Value Cutoff für Suchen |
| `min_identity` | 0.3 | Minimale Sequenzidentität |
| `max_intergenic_distance` | 50000 | Max. Abstand für Clustering |

---

## Output-Dateien

```
results/
├── report.html              # Hauptreport mit allen Ergebnissen
├── synteny_plot.html        # Interaktiver Syntenie-Plot
├── synteny_plot.png         # Statisches Bild
├── phylogenetic_tree.nwk    # Newick-Format Baum
├── all_proteins.faa         # Alle gefundenen Proteine
├── all_cds.fna              # Alle CDS-Sequenzen
├── alignment.faa            # MAFFT-Alignment
├── hits_summary.tsv         # Tabellarische Zusammenfassung
└── per_genome/
    ├── {genome_id}_hits.bed
    ├── {genome_id}_proteins.faa
    └── {genome_id}_annotations.gff
```

---

## Technische Abhängigkeiten

- **Nextflow** ≥ 24.x (DSL2)
- **MMseqs2**: Schnelle Sequenzsuche
- **Miniprot**: Protein-zu-Genom Alignment
- **Prodigal**: Ab-initio Genvorhersage
- **MAFFT**: Multiple Sequence Alignment
- **FastTree**: Phylogenetische Baumberechnung
- **Python 3.8+**: Für alle Skripte (ohne BioPython)
