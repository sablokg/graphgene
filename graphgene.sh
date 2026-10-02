#!/usr/bin/env bash
#
# ============================================================================
#  pangene_pipeline.sh Gaurav Sablok gsablok@proton.me
#  Complete pangenome GENE graph pipeline using pangene (Heng Li, 2024)
#  https://github.com/lh3/pangene   |   paper: Bioinformatics 2024, btae456
#
#  Unlike sequence-graph tools (pggb, minigraph-cactus, vg), pangene builds a
#  graph of GENES: nodes = marker genes, edges = genomic adjacency, walks =
#  genomes. It's built for exploring gene order/copy-number/presence-absence
#  across many assemblies (works for both eukaryotes and bacteria).
#
#  Pipeline stages:
#    1. Install pangene + miniprot + k8 (JS engine) + gfatools (viewer)
#    2. Build the input protein set (one canonical protein per gene)
#    3. Align that protein set to EVERY input genome with miniprot -> PAF
#    4. Build the gene graph from all PAFs                     -> pangene
#    5. Analyze: bubbles (gene-level variation) + presence/absence matrix
#    6. Optional: launch gfa-server for interactive visualization
#
#  Usage:
#    ./pangene_pipeline.sh -a annotation.gtf -r reference.fa \
#         -g genomes_dir/ -o results/
#
#    genomes_dir/ should contain one *.fna or *.fa per assembly
#    (the reference genome used for the annotation can be included too)
# ============================================================================

set -euo pipefail

# ---------------------------- CONFIG ----------------------------------------
GTF=""              # gene annotation (GTF/GFF3) for the protein set, e.g. GENCODE
REF_GENOME=""        # genome fasta matching $GTF (used to extract protein seqs)
GENOME_DIR=""        # directory of assemblies (*.fa/*.fna/*.fasta) to compare
OUTDIR="pangene_results"
THREADS=16
CANONICAL_ONLY=1     # 1 = one isoform per gene (-c), recommended for cleaner graphs
BACTERIA=0           # 1 = add -S to miniprot (disable splicing) for bacterial genomes

usage() {
  cat <<EOF
Usage: $0 -a annotation.gtf -r reference_genome.fa -g genomes_dir/ [-o outdir] [-t threads] [-b]
  -a  Gene annotation (GTF/GFF3), used to build the protein set (skip with -p)
  -r  Reference genome fasta matching -a
  -p  Provide a pre-built protein FASTA instead of -a/-r
  -g  Directory containing genome assemblies to compare (*.fa/*.fna)
  -o  Output directory (default: pangene_results)
  -t  Threads (default: 16)
  -b  Bacterial mode (disable splicing in miniprot)
EOF
  exit 1
}

PROTEIN_FAA=""
while getopts "a:r:p:g:o:t:bh" opt; do
  case $opt in
    a) GTF="$OPTARG" ;;
    r) REF_GENOME="$OPTARG" ;;
    p) PROTEIN_FAA="$OPTARG" ;;
    g) GENOME_DIR="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    b) BACTERIA=1 ;;
    h) usage ;;
    *) usage ;;
  esac
done

[[ -z "$GENOME_DIR" ]] && usage
[[ -z "$PROTEIN_FAA" && ( -z "$GTF" || -z "$REF_GENOME" ) ]] && usage

mkdir -p "$OUTDIR"/{bin,proteins,paf,graph,logs}
cd "$OUTDIR"

# ---------------------------- 1. INSTALL TOOLS -------------------------------
echo "=== [1/6] Installing pangene, miniprot, k8, gfatools ==="

install_tools() {
  cd bin

  # miniprot: protein-to-genome aligner (only aligner pangene supports)
  if [[ ! -x miniprot/miniprot ]]; then
    git clone https://github.com/lh3/miniprot 2>>../logs/install.log
    (cd miniprot && make -j"$THREADS") >>../logs/install.log 2>&1
  fi

  # pangene: the graph builder itself
  if [[ ! -x pangene/pangene ]]; then
    git clone https://github.com/lh3/pangene 2>>../logs/install.log
    (cd pangene && make -j"$THREADS") >>../logs/install.log 2>&1
  fi

  # k8: JS shell used to run pangene.js helper scripts (getaa, call, gfa2matrix)
  if [[ ! -x k8 ]]; then
    curl -L https://github.com/attractivechaos/k8/releases/download/v1.2/k8-1.2.tar.bz2 \
      | tar jxf - >>../logs/install.log 2>&1
    cp k8-1.2/k8-"$(uname | tr 'A-Z' 'a-z')"-x86_64 k8 2>/dev/null || \
      cp k8-1.2/k8-Linux k8  # fallback naming across k8 releases
    chmod +x k8
  fi

  # gfatools: provides gfa-server for interactive graph browsing (optional)
  if [[ ! -x gfatools/gfatools ]]; then
    git clone https://github.com/lh3/gfatools 2>>../logs/install.log
    (cd gfatools && make -j"$THREADS") >>../logs/install.log 2>&1
  fi

  cd ..
}
install_tools

MINIPROT=bin/miniprot/miniprot
PANGENE=bin/pangene/pangene
PANGENE_JS=bin/pangene/pangene.js
K8=bin/k8

export PATH="$PWD/bin:$PATH"

# ---------------------------- 2. BUILD PROTEIN SET ---------------------------
echo "=== [2/6] Building protein set ==="

if [[ -n "$PROTEIN_FAA" ]]; then
  cp "$PROTEIN_FAA" proteins/proteins.faa
else
  # getaa extracts one protein per transcript, named GENE:PROTEIN_ID
  # -c restricts to the canonical isoform per gene -> cleaner graphs
  if [[ "$CANONICAL_ONLY" -eq 1 ]]; then
    "$K8" "$PANGENE_JS" getaa -c "$GTF" "$REF_GENOME" > proteins/proteins.faa \
      2>logs/getaa.log
  else
    "$K8" "$PANGENE_JS" getaa "$GTF" "$REF_GENOME" > proteins/proteins.faa \
      2>logs/getaa.log
  fi
fi
N_PROT=$(grep -c '^>' proteins/proteins.faa)
echo "  -> $N_PROT protein sequences ready (proteins/proteins.faa)"

# ---------------------------- 3. ALIGN PROTEINS TO EACH GENOME ---------------
echo "=== [3/6] Aligning protein set to each genome (miniprot) ==="

MP_OPTS="--outs=0.97 --no-cs -Iut${THREADS}"
[[ "$BACTERIA" -eq 1 ]] && MP_OPTS="$MP_OPTS -S"

shopt -s nullglob
GENOME_FILES=("$GENOME_DIR"/*.fa "$GENOME_DIR"/*.fna "$GENOME_DIR"/*.fasta)
[[ ${#GENOME_FILES[@]} -eq 0 ]] && { echo "No genome files found in $GENOME_DIR"; exit 1; }

for genome in "${GENOME_FILES[@]}"; do
  name=$(basename "$genome" | sed -E 's/\.(fa|fna|fasta)(\.gz)?$//')
  paf="paf/${name}.paf"
  if [[ ! -s "$paf" ]]; then
    echo "  -> aligning $name"
    "$MINIPROT" $MP_OPTS "$genome" proteins/proteins.faa > "$paf" 2>>logs/miniprot.log
  fi
done

# ---------------------------- 4. BUILD THE PANGENE GRAPH ---------------------
echo "=== [4/6] Constructing pangene graph ==="

"$PANGENE" paf/*.paf > graph/pangene.gfa 2>logs/pangene.log
echo "  -> graph/pangene.gfa  ($(wc -l < graph/pangene.gfa) lines)"

# ---------------------------- 5. ANALYZE THE GRAPH ---------------------------
echo "=== [5/6] Analyzing graph: bubbles + presence/absence matrix ==="

# "Bibubbles": local subgraphs capturing gene order/copy-number/presence changes
"$K8" "$PANGENE_JS" call graph/pangene.gfa > graph/bubbles.txt 2>logs/call.log

# Binary presence/absence matrix (genes x genomes), Rtab format -> ready for
# downstream pangenome stats (e.g. core/accessory gene counts, PCA in R)
"$K8" "$PANGENE_JS" gfa2matrix graph/pangene.gfa > graph/gene_presence_absence.Rtab \
  2>logs/gfa2matrix.log

N_GENOMES=${#GENOME_FILES[@]}
N_CORE=$(awk -v n="$N_GENOMES" 'NR>1{c=0; for(i=2;i<=NF;i++) c+=($i>0); if(c==n) core++} END{print core+0}' \
  graph/gene_presence_absence.Rtab)
echo "  -> ${N_GENOMES} genomes compared; ${N_CORE} core genes (present in all)"

# ---------------------------- 6. VISUALIZATION (optional) --------------------
cat <<EOF > logs/README_visualize.txt
To interactively browse the graph:

  cd $OUTDIR
  bin/gfatools/gfatools ... (or use BandageNG on graph/pangene.gfa directly)

  # Or run pangene's own gfa-server (from a pangene release tarball, not built
  # by 'make'):
  curl -L https://github.com/lh3/pangene/releases/download/v1.1/pangene-1.1-bin.tar.bz2 | tar jxf -
  cd pangene-1.1-bin
  bin_linux-x64/gfa-server -d html ../graph/pangene.gfa.gz
  # open http://127.0.0.1:8000 and search for a gene name
EOF

echo "=== Done. Key outputs: ==="
echo "  graph/pangene.gfa                 - master gene graph (GFA)"
echo "  graph/bubbles.txt                 - gene-level structural variation calls"
echo "  graph/gene_presence_absence.Rtab   - presence/absence matrix for downstream stats"
echo "  logs/README_visualize.txt          - how to view the graph interactively"
