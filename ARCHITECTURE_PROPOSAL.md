# SynTerra Architecture Improvements - Vorschlag

## Aktuelle Probleme

### 1. Flanking Genes Qualität
- **Mit GFF**: Echte Gen-Annotationen (z.B. `gene-Melt`, `gene-LOC726827`)
- **Ohne GFF**: Prodigal-ORF-Predictions (`pred_OV788322.1_15592307`) → oft fragmentiert/ungenau

### 2. PHYLO_SORT Position
- **Aktuell**: Am Ende der Pipeline (sortiert nur Results)
- **Problem**: Genome werden unsortiert durchsucht → nahe verwandte Genome zuerst wäre besser für iterative Erweiterung

### 3. ANNOTATE_STRUCTURE Timing
- **Aktuell**: Einmalig nach ITERATIVE_SEARCH
- **Problem**: Gefundene Gene werden nicht sofort für weitere Suche genutzt

---

## Vorgeschlagene Architektur

```
                    ┌─────────────────┐
                    │   INPUT         │
                    │  - GOI Query    │
                    │  - Home Genom   │
                    │  - Target Taxa  │
                    └────────┬────────┘
                             │
              ┌──────────────┴──────────────┐
              ▼                             ▼
    ┌─────────────────┐           ┌─────────────────┐
    │  FETCH_QUERY    │           │ FETCH_RELATED   │
    │                 │           │    GENOMES      │
    └────────┬────────┘           └────────┬────────┘
             │                             │
             │                             ▼
             │                   ┌─────────────────┐
             │                   │  PHYLO_SORT     │  ◄── NEU: Am Anfang!
             │                   │ (MASH Distance) │
             │                   │ Sortiere nach   │
             │                   │ Verwandtschaft  │
             │                   └────────┬────────┘
             │                             │
             ▼                             │
    ┌─────────────────┐                    │
    │  LOCATE_GENE    │                    │
    └────────┬────────┘                    │
             │                             │
             ▼                             │
    ┌─────────────────────────────────────────────────┐
    │           INTELLIGENT_ANNOTATION                │
    │                                                 │
    │  Input: Home Genome + GOI Position              │
    │                                                 │
    │  1. Wenn GFF vorhanden:                         │
    │     → Extrahiere annotierte Flanking Genes      │
    │                                                 │
    │  2. Wenn KEINE GFF:                             │
    │     → Führe Miniprot mit Home-Proteom aus       │
    │     → Miniprot findet Exon-Boundaries           │
    │     → Supplement mit Prodigal für Lücken        │
    │                                                 │
    │  Output: Flanking Genes (hochwertig annotiert)  │
    └────────────────────────┬────────────────────────┘
                             │
                             ▼
    ╔═══════════════════════════════════════════════════════════════╗
    ║                   ITERATIVE_SEARCH (Wavefront)                ║
    ║                                                               ║
    ║   Für jedes Genom (sortiert nach phylogenetischer Distanz):   ║
    ║                                                               ║
    ║   ┌─────────────────────────────────────────────────────┐     ║
    ║   │  1. SEARCH: MMseqs2 tblastn-style                   │     ║
    ║   │     Query: Current expanded DB                       │     ║
    ║   │     Target: Genom                                    │     ║
    ║   └────────────────────┬────────────────────────────────┘     ║
    ║                        ▼                                      ║
    ║   ┌─────────────────────────────────────────────────────┐     ║
    ║   │  2. ANNOTATE: Miniprot für Treffer-Regionen         │     ║
    ║   │     - Findet exakte Exon-Boundaries                  │     ║
    ║   │     - Generiert splice-aware CDS                     │     ║
    ║   │     - Speichert GFF für spätere Referenz             │     ║
    ║   └────────────────────┬────────────────────────────────┘     ║
    ║                        ▼                                      ║
    ║   ┌─────────────────────────────────────────────────────┐     ║
    ║   │  3. VALIDATE: RBH Check + Domain Annotation         │     ║
    ║   │     - InterProScan/Pfam für Domänen                  │     ║
    ║   │     - Nutze Domänen als zusätzlichen Filter          │     ║
    ║   └────────────────────┬────────────────────────────────┘     ║
    ║                        ▼                                      ║
    ║   ┌─────────────────────────────────────────────────────┐     ║
    ║   │  4. EXPAND: Füge validierte Gene zur DB hinzu       │     ║
    ║   │     - Neue Gene = neue Queries für nächstes Genom   │     ║
    ║   └─────────────────────────────────────────────────────┘     ║
    ║                                                               ║
    ║   ───────────── Wiederhole für nächstes Genom ─────────────   ║
    ╚═══════════════════════════════════════════════════════════════╝
                             │
                             ▼
              ┌──────────────┴──────────────┐
              ▼                             ▼
    ┌─────────────────┐           ┌─────────────────┐
    │  COMPUTE_TREE   │           │  PLOT_SYNTENY   │
    │  (Final Tree)   │           │                 │
    └────────┬────────┘           └────────┬────────┘
             │                             │
             └──────────────┬──────────────┘
                            ▼
                  ┌─────────────────┐
                  │ GENERATE_REPORT │
                  └─────────────────┘
```

---

## Detaillierte Änderungen

### 1. PHYLO_SORT am Anfang (via MASH)

```python
# Neue Logik: Berechne MASH-Distanzen vor Suche
def calculate_genome_order(home_genome, target_genomes):
    """
    Verwendet MASH für schnelle Genom-Distanzberechnung.
    Sortiert Genome nach Distanz zum Home-Genom.
    """
    # MASH sketch für alle Genome
    # MASH dist home vs. alle targets
    # Sortiere: nächste Verwandte zuerst
```

**Vorteil**: Nahe verwandte Genome werden zuerst durchsucht → höhere Hit-Rate → bessere Query-Expansion für entfernte Genome.

### 2. Intelligente Annotation (Miniprot + Prodigal Hybrid)

```python
def intelligent_annotation(genome, goi_position, reference_proteome=None):
    """
    Kombiniert Miniprot und Prodigal für optimale Annotation.
    """
    
    # Schritt 1: Wenn Referenz-Proteom vorhanden
    if reference_proteome:
        # Miniprot: Mappt bekannte Proteine auf Genom
        miniprot_hits = run_miniprot(genome, reference_proteome)
        
        # Konvertiere zu GFF-ähnlicher Annotation
        for hit in miniprot_hits:
            # Exakte Exon-Boundaries von Miniprot
            gene = {
                'id': hit['target_id'],
                'exons': hit['cds_parts'],
                'strand': hit['strand'],
                'source': 'miniprot'
            }
    
    # Schritt 2: Prodigal für Lücken
    prodigal_orfs = run_prodigal(genome)
    
    # Schritt 3: Merge (Miniprot hat Priorität)
    final_genes = merge_annotations(miniprot_hits, prodigal_orfs)
    
    return final_genes
```

### 3. Integrierte Domänen-Annotation in ITERATIVE_SEARCH

```python
def annotate_and_validate(protein_sequence, domains_db):
    """
    Annotiert Domänen und nutzt sie als zusätzlichen Validierungsfilter.
    """
    
    # hmmscan gegen Pfam
    domains = run_hmmscan(protein_sequence, domains_db)
    
    # Erwartete Domänen für Gen-Familie prüfen
    expected_domains = get_expected_domains(gene_family)
    
    if has_expected_domains(domains, expected_domains):
        return True, domains
    else:
        return False, domains
```

---

## Vorteile der neuen Architektur

1. **Bessere Sensitivität**: 
   - Phylogenetisch sortierte Suche startet mit nahen Verwandten
   - Gefundene Orthologe erweitern Query-Pool für entfernte Genome

2. **Höhere Präzision**:
   - Miniprot für exakte Exon-Boundaries statt Prodigal-ORFs
   - Domänen-basierte Validierung als zusätzlicher Filter

3. **Bessere Annotation**:
   - Jedes Genom bekommt sofort hochwertige Gen-Annotation
   - GFF-Dateien für alle Treffer

4. **Schnellere Konvergenz**:
   - Frühe Treffer verbessern Suche in späteren Genomen

---

## Implementierungs-Prioritäten

### Phase 1: Phylogenetische Sortierung
- [ ] MASH-basierte Distanzberechnung vor ITERATIVE_SEARCH
- [ ] Genome nach Distanz sortieren

### Phase 2: Verbesserte Home-Annotation
- [ ] Miniprot-basierte Flanking Gene Annotation wenn kein GFF
- [ ] Hybrid mit Prodigal für Lücken

### Phase 3: Integrierte Domänen-Annotation
- [ ] hmmscan in ITERATIVE_SEARCH Loop
- [ ] Domänen-basierte Filterung
