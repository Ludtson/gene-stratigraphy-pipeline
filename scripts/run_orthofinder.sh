#!/usr/bin/env bash
set -euo pipefail

# ==========================================
# Orthology & Lineage Analysis Pipeline
# Version: 1.0.0
# Author: Adekola Owoyemi (Casola Lab, Texas A&M University)
# ==========================================

# --- SCRIPT DESCRIPTION ---
#
# run_orthofinder.sh
#
# This script runs OrthoFinder on a collection of protein FASTA files
# (typically produced from Step 1 of the pipeline). It prepares input
# data, optionally renames species using a mapping file, and executes
# orthology inference.
#
# Workflow:
# 1. Validate input data and environment
# 2. Prepare canonical FASTA directory (protein/)
# 3. Create temporary working directory
# 4. Symlink FASTA files for OrthoFinder execution
# 5. Run OrthoFinder (DIAMOND-based)
# 6. Move results into structured output directory
#
# Outputs:
#   output_dir/
#     ├── protein/         # canonical FASTA files (persistent)
#     ├── orthofinder/     # OrthoFinder results
#     └── logs/            # execution logs
#
# Notes:
# - Conda environment must be activated before running this script.
# - Map file is optional. If not provided, species names are derived from
#   FASTA filenames, stripping a trailing "_final" suffix if present (Step
#   1's fixed output naming convention: <species>_final.faa under
#   5_final_proteins/). This keeps inferred species names consistent with
#   a rooted species tree built from clean species names.
#
# Usage:
#   bash run_orthofinder.sh -i <input_dir> -o <output_dir> [-t threads] [-m map_file]
#
# Example:
#   bash run_orthofinder.sh -i step1_output -o step2_output -t 16
#


# ==============================================================================
# DEFAULTS
# ==============================================================================
threads=8
map_file=""
log_file="orthofinder.log"
species_tree=""


# ==============================================================================
# LOGGING
# ==============================================================================
log() {
    local level="$1"
    local msg="$2"
    echo "$(date '+%Y-%m-%d %H:%M:%S') [$level] $msg"
}


# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -i|--input) input_dir="$2"; shift ;;
        -o|--output) output_dir="$2"; shift ;;
        -t|--threads) threads="$2"; shift ;;
        -m|--map) map_file="$2"; shift ;;
        -s|--species-tree) species_tree="$2"; shift ;;
        -h|--help)
            echo "
Usage:
  $0 -i <input_dir> -o <output_dir> [-t threads] [-m map_file]

Description:
  Run OrthoFinder on a set of protein FASTA files.

Required:
  -i, --input     Directory containing protein FASTA files (Step 1 output)
  -o, --output    Output directory for results

Optional:
  -t, --threads   Number of threads (default: 8)
  -m, --map       Tab-separated file for renaming species
  -s, --species-tree  Rooted Newick species tree file

Map file format (optional):
  Species    Basename
  Homo_sapiens    hsap
  Mus_musculus    mmus

  'Basename' should match the prefix of FASTA files in input_dir.

Output:
  output_dir/
    ├── protein/        # prepared FASTA files
    ├── orthofinder/    # OrthoFinder results
    └── logs/           # log files

Example:
  $0 -i step1_output -o step2_output -t 16
"
            exit 0
            ;;
        *)
            log ERROR "Unknown parameter: $1"
            exit 1
            ;;
    esac
    shift
done


# ==============================================================================
# VALIDATION
# ==============================================================================
[[ -z "${input_dir:-}" || -z "${output_dir:-}" ]] && {
    log ERROR "Missing required arguments. Use -h for help."
    exit 1
}

[[ ! -d "$input_dir" ]] && {
    log ERROR "Input directory not found: $input_dir"
    exit 1
}

if [[ -n "$species_tree" ]]; then
    if [[ ! -f "$species_tree" ]]; then
        log ERROR "Species tree file not found: $species_tree"
        exit 1
    fi
    species_tree=$(realpath "$species_tree")
fi

command -v orthofinder >/dev/null 2>&1 || {
    log ERROR "orthofinder not found in PATH"
    exit 1
}

# Resolve absolute paths
input_dir=$(realpath "$input_dir")
output_dir=$(realpath "$output_dir")
[[ -n "$map_file" ]] && map_file=$(realpath "$map_file")


# ==============================================================================
# SETUP DIRECTORIES
# ==============================================================================
mkdir -p "$output_dir/protein"
mkdir -p "$output_dir/orthofinder"
mkdir -p "$output_dir/logs"

log_file="$output_dir/logs/$log_file"

# Clean previous protein files (safe reruns)
rm -f "$output_dir/protein"/*.fa* 2>/dev/null || true


# ==============================================================================
# PREPARE FASTA FILES
# ==============================================================================
log INFO "Preparing FASTA files..."

shopt -s nullglob

if [[ -n "$map_file" ]]; then
    log INFO "Using species map file: $map_file"

    tail -n +2 "$map_file" | tr -d '\r' | while IFS=$'\t' read -r species basename; do
        [[ -z "$species" || -z "$basename" ]] && continue

        # Robust file matching (handles extension variation)
        src=$(ls "$input_dir"/"$basename"*.fa* 2>/dev/null | head -n 1)

        if [[ -z "$src" ]]; then
            log ERROR "Missing FASTA for basename: $basename"
            exit 1
        fi

        dest="$output_dir/protein/${species}.faa"
        cp "$src" "$dest"
    done

else
    log WARN "No map file provided — deriving species names from filenames"

    for f in "$input_dir"/*.fa "$input_dir"/*.faa "$input_dir"/*.fasta; do
        fname=$(basename "$f")
        ext="${fname##*.}"
        base="${fname%.*}"
        # Step 1 (protein-preprocessing-isoform-pipeline) writes
        # <species>_final.<ext> under 5_final_proteins/ -- strip that
        # suffix so the inferred species name matches a rooted species
        # tree built from clean species names (e.g. via -s).
        species="${base%_final}"
        cp "$f" "$output_dir/protein/${species}.${ext}"
    done
fi

log INFO "FASTA preparation complete"


# ==============================================================================
# TEMP DIRECTORY
# ==============================================================================
temp_dir=$(mktemp -d -t orthofinder_XXXX)
trap 'rm -rf "$temp_dir"' EXIT

shopt -s nullglob

files=(
    "$output_dir"/protein/*.fa
    "$output_dir"/protein/*.faa
    "$output_dir"/protein/*.fasta
    "$output_dir"/protein/*.fas
)

if [ ${#files[@]} -eq 0 ]; then
    log ERROR "No FASTA files found in $output_dir/protein"
    exit 1
fi

ln -s "${files[@]}" "$temp_dir"/


# ==============================================================================
# RUN ORTHOFINDER
# ==============================================================================
log INFO "Running OrthoFinder with $threads threads..."

if [[ -n "$species_tree" ]]; then
    log INFO "Using provided species tree: $species_tree"
    
    orthofinder -f "$temp_dir" -t "$threads" -S diamond -s "$species_tree" \
    >> "$log_file" 2>&1
else
    log INFO "No species tree provided — inferring tree"
    
    orthofinder -f "$temp_dir" -t "$threads" -S diamond \
    >> "$log_file" 2>&1
fi


# ==============================================================================
# MOVE RESULTS
# ==============================================================================
log INFO "Collecting results..."

if [[ -d "$temp_dir/OrthoFinder" ]]; then
    cp -r "$temp_dir/OrthoFinder/"* "$output_dir/orthofinder/"
else
    log ERROR "OrthoFinder output not found"
    exit 1
fi


# ==============================================================================
# DONE
# ==============================================================================
log INFO "OrthoFinder run completed successfully"
log INFO "Results directory: $output_dir/orthofinder"
log INFO "Log file: $log_file"
