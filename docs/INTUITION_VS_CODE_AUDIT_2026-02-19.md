# SynTerra Audit: Intuition vs. tatsaechlicher Code (2026-02-19)

## Ziel
Diese Analyse prueft deine Soll-Logik Schritt fuer Schritt gegen den aktuellen Codepfad und gegen die erzeugten Snake-Outputs (`results_3snake_3ftx_v7_docker`).

Fokus: nur Punkte, die biologisch wirklich relevant sind (TOP1MT/Ly6/3FTX-Kontext, Annotationstreue, Locus-Logik, GOI-Erkennung).

## Protokoll (was ich konkret geprueft habe)
1. Workflow-Orchestrierung in `main.nf` (Easy/Pro-Mode, Prozessreihenfolge, Locus-Split, Flanking, iterative Suche, Clustering, Plot).
2. Kernmodule/Skripte fuer die Soll-Schritte:
   - Input-Aufloesung: `bin/resolve_gene_input.py`, `modules/resolve_query.nf`
   - Home/Target-Genome-Auswahl: `bin/fetch_home_genome.py`, `bin/fetch_related_genomes.py`
   - GOI-Lokalisierung/Locus-Split: `modules/locate_gene.nf`, `bin/merge_hits.py`, `bin/split_loci.py`
   - Flanking-Extraktion/Annotation: `bin/extract_flanking_genes.py`
   - GOI-Exonannotation: `bin/annotate_goi_exons.py`
   - Iterative Suche/Region-Modellierung: `bin/iterative_search_runner.py`, `bin/cluster_grs.py`
   - Plot/Labeling: `bin/plot_synteny.py`
3. Ergebnisdateien gegengeprueft:
   - `results_3snake_3ftx_v7_docker/intermediate/annotate_goi/goi_info.json`
   - `results_3snake_3ftx_v7_docker/intermediate/annotate_goi/goi_exons.faa`
   - `results_3snake_3ftx_v7_docker/intermediate/flanking/synteny_block_locus_1.bed`
   - `results_3snake_3ftx_v7_docker/home_genome/home_genome/home_genome.gff`
   - `results_3snake_3ftx_v7_docker/downloaded_genomes/easy_mode_genomes/*`

## Soll vs Ist (deine Punkte)

### 1) Input: UniProt/NCBI/rohe Sequenz + Home-Species
**Soll:** flexible Eingabe, Home-Spezies ableitbar oder explizit.

**Ist:** weitgehend erfuellt.
- Easy mode unterstuetzt lokale FASTA direkt oder Aufloesung via Resolver (`main.nf:162-183`).
- Resolver erkennt Datei/UniProt/NCBI/Symbol (`bin/resolve_gene_input.py:34-53`, `bin/resolve_gene_input.py:313-330`).
- Spezies wird uebernommen/abgeleitet, sonst harter Fehler (`main.nf:185-193`).

**Bewertung:** passt im Kern.

### 2) GOI im Home-Genom lokalisieren, inkl. Duplikationen
**Soll:** GOI robust finden, Mehrfachloci trennen.

**Ist:** teilweise erfuellt, aber zu permissiv.
- Lokalisierung via MMseqs + BLAST in `LOCATE_GENE` (`modules/locate_gene.nf:17-77`).
- Hits werden mit 1 kb Gap zusammengefuehrt (`bin/merge_hits.py:78-89`).
- Locus-Split trennt nach Distanz (`bin/split_loci.py:58-68`).

**Abweichung:**
- `merge_hits.py` filtert standardmaessig mit `--max_evalue 1e-3` (`bin/merge_hits.py:11-13`), aber `LOCATE_GENE` uebergibt keinen Wert (`modules/locate_gene.nf:77`).
- Das ist lockerer als die eigentliche Suche und kann bei kurzen Toxinen mehr Rauschen in Loci bringen.

### 3) Genome-Auswahl (easy mode)
**Soll:** gute Auswahl-Logik; vorgegebene Targets beruecksichtigen.

**Ist:** grundsaetzlich gut.
- Ranking mit RefSeq/Assembly-Level/Contigs/N50 etc. (`bin/fetch_related_genomes.py:311-361`).
- Bad-quality policy vorhanden (`bin/fetch_related_genomes.py:375-395`).
- Download versucht GFF + RefSeq-Fallback (`bin/fetch_related_genomes.py:668-739`).

**Abweichung:**
- GFF ist nicht garantiert. In deinem Run hat ein Target nur `.fna` (kein `.gff`):
  - vorhanden: `GCA_023653725.1.fna`
  - fehlend: `GCA_023653725.1.gff`
  - siehe `results_3snake_3ftx_v7_docker/downloaded_genomes/easy_mode_genomes`.

### 4) Flanking-Gene an Loci (nicht pro Einzelhit redundant)
**Soll:** bei Tandem-Regionen Locus als Ganzes flankieren, nicht redundant pro Hit.

**Ist:** abweichend.
- `extract_flanking_genes.py` waehlt pro Zielregion den naechsten Gen-Zentrumspunkt und nimmt `n_flank` drum herum (`bin/extract_flanking_genes.py:339-364`).
- Bei mehreren Regionen im selben Locus wird mehrfach Fenster-Logik angewendet; danach nur ID-dedupliziert (`bin/extract_flanking_genes.py:397-405`).

**Konsequenz:**
- Kein explizites "Locus-Grenzen zuerst, dann einmal flankieren"-Verhalten.
- Redundanz/Verzerrung moeglich, besonders bei Cluster von GOI-Duplikationen.

### 5) Annotation verwenden (CDS/Exons) + gene names statt kryptischer IDs
**Soll:** CDS/Exon-basiert und fuer Menschen lesbare Gen-Namen.

**Ist:** halb erfuellt.
- Mit GFF werden CDS/Exons wirklich genutzt (`bin/extract_flanking_genes.py:419-478`).
- Aber IDs/Labels kommen primaer aus `attrs['ID']` (`bin/extract_flanking_genes.py:400`), nicht aus `gene`/`Name` als primaerer Ausgabe-ID.
- Ergebnis ist entsprechend kryptisch: `gene-E2320_...` in `synteny_block_locus_1.bed`.

**Konsequenz:**
- Plot wirkt biologisch "blind", weil z. B. `TOP1MT` visuell nicht auftaucht, selbst wenn semantisch nahe Features da waeren.

### 6) Wave-Search: erst Flanking, dann GOI-sensitiv im Locus
**Soll:** Flanking verankert Regionen, danach sensitiver GOI-Search (inkl. SW/tblastn) lokal begrenzt.

**Ist:** konzeptionell weitgehend passend.
- Block-Seeding bevorzugt flanking hits (`bin/iterative_search_runner.py:2503-2510`).
- Danach GOI-fokussierte augmentierte Suche im gepaddeten Block (`bin/iterative_search_runner.py:1726-1733`, `bin/iterative_search_runner.py:1648-1654`).
- SW/tblastn-Rescue ist integriert.

**Kritischer Bruch (siehe unten):** GOI-Exon-IDs koennen falsch klassifiziert werden.

## Schwere Fehler (priorisiert)

## CRITICAL-1: GOI-Exons verlieren GOI-Praefix im Miniprot/Hit-Modus
- In hit-basierter GOI-Annotation werden Exons als `exon_1`, `exon_2`, ... erzeugt (`bin/annotate_goi_exons.py:1030`).
- GOI-Klassifikation im Iterationskern erkennt GOI aber ueber Praefix `GOI_` (`bin/iterative_search_runner.py:228-242`).
- Damit fallen diese Exons als "nicht GOI" durch.

**Run-Beleg:**
- `results_3snake_3ftx_v7_docker/intermediate/annotate_goi/goi_info.json` zeigt `"id": "exon_1" ...`.
- `results_3snake_3ftx_v7_docker/intermediate/annotate_goi/goi_exons.faa` enthaelt `>GOI_P60615` plus `>exon_1`, `>exon_2`, `>exon_3`.

**Biologische Folge:**
- GOI-Signal wird in Folgeprozessen teilweise wie Flanking/sonstige Query behandelt.
- Genauigkeit der GOI-Verfolgung ueber Genome sinkt.

## CRITICAL-2: Flanking-Locus-Logik nicht locus-zentriert
- Flanking-Fenster wird pro Region um den jeweils "naechsten Gene-Center" gezogen (`bin/extract_flanking_genes.py:343-364`).
- Gewuenscht waere: pro Locus erst Gesamtspanne definieren, dann ein gemeinsamer Flanking-Block.

**Biologische Folge:**
- Redundanz und potenziell inkonsistente Nachbarschaftsdefinition bei Duplikationsclustern.
- Schwaecht direkte Vergleichbarkeit mit Ground Truth (z. B. TOP1MT-Ly6-Kontext).

## HIGH-1: Gene-Labels sind technisch, nicht biologisch lesbar
- Ausgabe verwendet ueberwiegend `ID` statt symbolischer Namen (`bin/extract_flanking_genes.py:400`).
- In Bed/Plot landen damit primaer `gene-E2320_...`-Labels.

**Biologische Folge:**
- Plot ist schwer interpretierbar und kaum paper-ready.

## HIGH-2: TOP1MT-Anker ist nirgends als harte Ziel-Bedingung verdrahtet
- Es gibt derzeit keine "TOP1MT als Flank-Anker"-Regel im Workflow.
- In deinem Home-GFF wurden 3FTX-Hits im Bereich `CM019150.1` gefunden (`home_genome.gff:191296-191305`),
  waehrend "putative DNA topoisomerase protein" an anderem Contig liegt (`home_genome.gff:368089ff`, `CM019165.1`).

**Biologische Folge:**
- Erwartete Story (TOP1MT -> Ly6 -> 3FTX) wird nicht gezielt getestet/erzwingbar gemacht.

## HIGH-3: Region-Scoring/"p-value" ist Platzhalter
- `cluster_grs.py` nutzt `p_value = 1.0 - score` (`bin/cluster_grs.py:362-364`, `bin/cluster_grs.py:431`).
- Das ist kein statistisch belastbarer Signifikanztest.

**Biologische Folge:**
- Confidence-Aussagen sind fuer Publikation nur eingeschraenkt belastbar.

## MEDIUM-1: E-Value-Logik inkonsistent zwischen Search und Merge
- Suche nutzt `params.search_evalue`, Merge standardmaessig `1e-3` ohne Param-Uebergabe (`modules/locate_gene.nf:41`, `modules/locate_gene.nf:77`, `bin/merge_hits.py:11-13`).

**Biologische Folge:**
- Niedrigqualitative Hits koennen beim Mergen unnoetig im Locus landen.

## Ground-Truth-Abgleich (TOP1MT/Ly6/3FTX)
1. 3FTX-annotierte Features sind im gewaehlten Home-Locus klar vorhanden (`CM019150.1`, z. B. `home_genome.gff:191297`, `home_genome.gff:191301`).
2. Explizites `TOP1MT` taucht in diesem Locus nicht auf; topoisomerase-Annotation liegt separat (`home_genome.gff:368089ff`, `CM019165.1`).
3. Flanking-Ausgabe ist primaer locus_tag/ID-basiert (`gene-E2320_...`), daher visuell kaum als TOP1MT/Ly6 interpretierbar.

## Fazit zur Kernfrage
Deine Intuition ist in der Gesamtarchitektur grundsaetzlich getroffen (Input -> GOI-Lokalisierung -> Flanking -> iterative Suche -> Region/Plot), aber an drei entscheidenden Stellen weicht die Implementierung stark ab:
1. GOI-Exon-ID-Handling (kritischer funktionaler Fehler),
2. Flanking nicht sauber locus-zentriert,
3. Annotation-/Namensdarstellung nicht biologiezentriert.

Damit ist der aktuelle Plotstand fuer die gewuenschte TOP1MT/Ly6/3FTX-Hypothese noch nicht robust/paper-ready.

## Priorisierte Fix-Strategie (Nature-ready, nur high-impact)
1. **GOI-ID-Konsistenz reparieren (CRITICAL-1).**
   - Alle GOI-Exons aus allen Pfaden konsistent als `GOI_<query>|exon_N` serialisieren.
   - Acceptance: keine `exon_N` ohne `GOI_` mehr in `goi_exons.faa`.
2. **Locus-zentrierte Flanking-Extraktion (CRITICAL-2).**
   - Erst Locusspanne aggregieren, dann einmal flankieren.
   - Acceptance: bei einem Tandem-Cluster genau ein konsistenter Flanking-Block.
3. **Biologische Label-Prioritaet in BED/Plot (HIGH-1).**
   - Label-Reihenfolge: `gene` > `Name` > `locus_tag` > `ID`.
   - Acceptance: TOP1MT/Ly6 (falls annotiert) erscheinen direkt im Plot.
4. **TOP1MT-Anker als optionaler Modus (HIGH-2).**
   - Optionaler Constraint/Report: "enthaltener Flankanker TOP1MT ja/nein".
5. **Score-Validitaet verbessern (HIGH-3).**
   - Platzhalter-p-value entfernen oder als "heuristic" klar kennzeichnen.

## Nicht-Ziele in dieser Runde
- Keine kosmetischen Refactors.
- Keine Aenderungen ohne klaren biologischen Mehrwert.
