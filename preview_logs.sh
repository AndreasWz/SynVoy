#!/bin/bash
# Quick preview of the new SynTerra pipeline logging

# ANSI color codes
c_reset="\033[0m"
c_dim="\033[2m"
c_red="\033[0;31m"
c_green="\033[0;32m"
c_yellow="\033[0;33m"
c_blue="\033[0;34m"
c_cyan="\033[0;36m"
c_white="\033[0;37m"

# ASCII Banner
echo -e "
${c_cyan}╔═══════════════════════════════════════════════════════════════╗
║                                                               ║
║   ${c_white}███████╗██╗   ██╗███╗   ██╗████████╗███████╗██████╗ ██████╗  ███████╗${c_cyan}  ║
║   ${c_white}██╔════╝╚██╗ ██╔╝████╗  ██║╚══██╔══╝██╔════╝██╔══██╗██╔══██╗██╔══██║${c_cyan}  ║
║   ${c_white}███████╗ ╚████╔╝ ██╔██╗ ██║   ██║   █████╗  ██████╔╝██████╔╝███████║${c_cyan}  ║
║   ${c_white}╚════██║  ╚██╔╝  ██║╚██╗██║   ██║   ██╔══╝  ██╔══██╗██╔══██╗██╔══██║${c_cyan}  ║
║   ${c_white}███████║   ██║   ██║ ╚████║   ██║   ███████╗██║  ██║██║  ██║██║  ██║${c_cyan}  ║
║   ${c_white}╚══════╝   ╚═╝   ╚═╝  ╚═══╝   ╚═╝   ╚══════╝╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝${c_cyan}  ║
║                                                               ║
║   ${c_white}Phylogenetically-informed syntenic ortholog discovery${c_cyan}     ║
║   ${c_dim}v2.0 | github.com/AndreasWz/SynTerra${c_cyan}                      ║
║                                                               ║
╚═══════════════════════════════════════════════════════════════╝${c_reset}
"

# Configuration
echo -e "${c_blue}═══════════════════════════════════════════════════════════════
${c_white}RUN CONFIGURATION${c_reset}
${c_blue}═══════════════════════════════════════════════════════════════${c_reset}
${c_dim}Query Gene      :${c_reset} ${c_green}tetramorium_query.fasta${c_reset}
${c_dim}Home Genome     :${c_reset} ${c_green}tetramorium_home.fasta${c_reset}
${c_dim}Mode            :${c_reset} ${c_yellow}standard${c_reset}
${c_dim}Flanking Genes  :${c_reset} 5
${c_dim}MMseqs Sens.    :${c_reset} 7.5
${c_dim}Output Dir      :${c_reset} results/test_run
${c_blue}═══════════════════════════════════════════════════════════════${c_reset}
"

echo -e "${c_cyan}⚡ Starting pipeline execution...${c_reset}\n"

# Phase 1
echo -e "${c_blue}═══════════════════════════════════════════════════════════════
${c_white}PHASE 1: Gene Location & Genome Preparation${c_reset}
${c_blue}═══════════════════════════════════════════════════════════════${c_reset}"

sleep 0.5
echo -e "${c_cyan}📥 Fetching related genomes...${c_reset}"
sleep 0.5
echo -e "${c_green}✓ Found 4 genomes from NCBI${c_reset}"
sleep 0.5
echo -e "${c_cyan}📂 Staging genome files...${c_reset}"
sleep 0.5
echo -e "${c_green}✓ Genomes staged successfully${c_reset}"
sleep 0.5
echo -e "${c_cyan}🔍 Locating query gene in home genome...${c_reset}"
sleep 0.5
echo -e "${c_green}✓ Gene located: chr1:245678-246890${c_reset}"
sleep 0.5
echo -e "${c_cyan}🧬 Extracting flanking genes...${c_reset}"
sleep 0.5
echo -e "${c_green}✓ Extracted 5 flanking genes${c_reset}"

# Phase 2
echo -e "\n${c_blue}═══════════════════════════════════════════════════════════════
${c_white}PHASE 2: Phylogenetic Ordering & Iterative Search${c_reset}
${c_blue}═══════════════════════════════════════════════════════════════${c_reset}"

sleep 0.5
echo -e "${c_cyan}🌳 Sorting genomes by phylogenetic distance...${c_reset}"
sleep 0.5
echo -e "${c_green}✓ Phylogenetic ordering complete for locus_1${c_reset}"
sleep 0.5
echo -e "${c_cyan}📊 Assessing genome quality...${c_reset}"
sleep 0.5
echo -e "${c_green}✓ Quality assessment complete${c_reset}"
sleep 0.5
echo -e "${c_cyan}🏠 Preparing home proteome database...${c_reset}"
sleep 0.5
echo -e "${c_green}✓ Home proteome database ready${c_reset}"
sleep 0.5
echo -e "${c_cyan}🔍 Running iterative phylogenetic search...${c_reset}"
sleep 1
echo -e "${c_green}✓ Iterative search complete: locus_1${c_reset}"

# Phase 3
echo -e "\n${c_blue}═══════════════════════════════════════════════════════════════
${c_white}PHASE 3: Region Clustering & Augmented Search${c_reset}
${c_blue}═══════════════════════════════════════════════════════════════${c_reset}"

sleep 0.5
echo -e "${c_cyan}🧩 Clustering genomic regions by synteny...${c_reset}"
sleep 0.7
echo -e "${c_green}✓ Clustered regions for locus_1 in genome_001${c_reset}"
sleep 0.3
echo -e "${c_green}✓ Clustered regions for locus_1 in genome_002${c_reset}"
sleep 0.3
echo -e "${c_green}✓ Clustered regions for locus_1 in genome_003${c_reset}"
sleep 0.5
echo -e "${c_cyan}🔬 Running augmented orthology search...${c_reset}"
sleep 0.7
echo -e "${c_green}✓ Augmented search complete for region genome_001_locus_1${c_reset}"
sleep 0.3
echo -e "${c_green}✓ Augmented search complete for region genome_002_locus_1${c_reset}"
sleep 0.3
echo -e "${c_green}✓ Augmented search complete for region genome_003_locus_1${c_reset}"

# Phase 4
echo -e "\n${c_blue}═══════════════════════════════════════════════════════════════
${c_white}PHASE 4: Phylogenetics & Visualization${c_reset}
${c_blue}═══════════════════════════════════════════════════════════════${c_reset}"

sleep 0.5
echo -e "${c_cyan}🌳 Computing phylogenetic trees...${c_reset}"
sleep 0.8
echo -e "${c_green}✓ Phylogenetic tree computed: locus_1${c_reset}"
sleep 0.5
echo -e "${c_cyan}📊 Generating synteny visualizations...${c_reset}"
sleep 0.8
echo -e "${c_green}✓ Synteny visualization complete${c_reset}"

# Phase 5
echo -e "\n${c_blue}═══════════════════════════════════════════════════════════════
${c_white}PHASE 5: Report Generation${c_reset}
${c_blue}═══════════════════════════════════════════════════════════════${c_reset}"

sleep 0.5
echo -e "${c_cyan}📝 Generating comprehensive report...${c_reset}"
sleep 0.8
echo -e "${c_green}✓ Analysis report generated successfully${c_reset}"

# Completion
echo -e "\n${c_blue}═══════════════════════════════════════════════════════════════${c_reset}"
echo -e "${c_green}✓ Pipeline completed successfully!${c_reset}"
echo -e "${c_blue}═══════════════════════════════════════════════════════════════${c_reset}"
echo -e "${c_white}  Results Directory:${c_reset} ${c_cyan}results/test_run${c_reset}"
echo -e "${c_white}  Duration:${c_reset}          ${c_dim}3m 42s${c_reset}"
echo -e "${c_white}  Tasks Completed:${c_reset}   ${c_dim}24${c_reset}"
echo -e "${c_blue}═══════════════════════════════════════════════════════════════${c_reset}"
echo -e "${c_dim}  Key outputs:${c_reset}"
echo -e "${c_dim}    • synterra_report.json    (Analysis summary)${c_reset}"
echo -e "${c_dim}    • synteny_plot.pdf        (Visualization)${c_reset}"
echo -e "${c_dim}    • expanded_databases/     (Ortholog databases)${c_reset}"
echo -e "${c_dim}    • augmented_regions/      (Identified regions)${c_reset}"
echo -e "${c_blue}═══════════════════════════════════════════════════════════════${c_reset}\n"
