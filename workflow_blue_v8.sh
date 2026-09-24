#!/usr/bin/env bash
set -Eeuo pipefail
export LC_ALL=C

# GWAS-submit BLUE workflow v8
#
# Single phenotype:
#   immediately delegates to the original workflow_v2.sh.
#
# >=2 phenotype files:
#   prepare:
#     - unpack files
#     - direct conversion to lme4 CSV
#     - run the one BLUE mode automatically determined by the environment structure
#     - preserve historical raw BLUE txt/csv in lme4/blue/
#     - combine normalized plan columns for one downstream GWAS preparation
#     - delegate to patched original workflow_v2.sh
#   admin/final:
#     directly delegate to original workflow_v2.sh.
#
# The patched original workflow is uploaded as workflow_v2.sh by server_blue.py.

HERE=$(cd "$(dirname "$0")" && pwd)
ORIGINAL="$HERE/workflow_v2.sh"

stage=${1:-}
pop=${2:-}
folder=${3:-}
filename=${4:-}
runid=${5:-}
tag=${6:-}

BASE=/data9/home/yzhao/GWAS_IRGSP1.0
staging="$BASE/.gwas-web/$runid"
bundle="$staging/phenotype.upload"
MAGIC='#GWAS_BLUE_BUNDLE_V4'

mkdir -p "$staging"
exec >> "$staging/run.log" 2>&1

finish_wrapper() {
    rc=$?
    trap - EXIT
    if [[ ! -f "$staging/$stage.exit" ]]; then
        printf '%s\n' "$rc" > "$staging/$stage.exit.tmp"
        mv "$staging/$stage.exit.tmp" "$staging/$stage.exit"
    fi
}
trap finish_wrapper EXIT

echo "============================================================"
echo "BLUE_WRAPPER_V8_START"
echo "STAGE: $stage"
echo "RUN_ID: $runid"
echo "POPULATION: $pop"
echo "PROJECT: $folder"
echo "============================================================"

[[ -s "$ORIGINAL" ]] || {
    echo "ERROR: missing remote workflow_v2.sh: $ORIGINAL"
    exit 1
}

# Single phenotype or later stages use upstream workflow.
if [[ $stage != prepare ]] || [[ ! -s $bundle ]] || \
   [[ $(head -n 1 "$bundle") != "$MAGIC" ]]; then
    echo "BLUE_WRAPPER_DELEGATE: original workflow_v2.sh"
    trap - EXIT
    exec bash "$ORIGINAL" "$@"
fi

[[ $pop =~ ^(GPall|GPallInd|GPallJap)$ ]] || {
    echo "ERROR: invalid population: $pop"
    exit 64
}

work="$BASE/$folder"
mkdir -p \
    "$work/phenotype" \
    "$work/lme4" \
    "$work/lme4/blue" \
    "$staging/blue-unpack"

echo "PROJECT_DIR_CREATED: $work"

# This marker tells the patched original workflow not to prepend the global
# tag "BLUE" to columns that already carry their own BLUE plan label.
printf 'BLUE\n' > "$staging/blue-mode.flag"

b64decode() {
    printf '%s' "$1" | base64 -d
}

declare -A file_name file_prefix file_loc file_year file_b64
declare -A plan_id plan_kind plan_label plan_output plan_title plan_indexes
traits=()
population_meta=""

while IFS=$'\t' read -r kind a b c d e f g; do
    case "$kind" in
        "$MAGIC")
            ;;
        META)
            value=$(b64decode "$b")
            [[ $a == population ]] && population_meta=$value
            ;;
        TRAIT)
            traits+=("$(b64decode "$b")")
            ;;
        FILE)
            idx=$a
            file_name[$idx]=$(b64decode "$b")
            file_prefix[$idx]=$(b64decode "$c")
            file_loc[$idx]=$(b64decode "$d")
            file_year[$idx]=$(b64decode "$e")
            file_b64[$idx]="$staging/blue-unpack/$idx.b64"
            : > "${file_b64[$idx]}"
            ;;
        DATA)
            printf '%s' "$b" >> "${file_b64[$a]}"
            ;;
        END)
            idx=$a
            base64 -d "${file_b64[$idx]}" \
                > "$staging/blue-unpack/$idx.raw"
            ;;
        PLAN)
            idx=$a
            plan_id[$idx]=$(b64decode "$b")
            plan_kind[$idx]=$(b64decode "$c")
            plan_label[$idx]=$(b64decode "$d")
            plan_output[$idx]=$(b64decode "$e")
            plan_title[$idx]=$(b64decode "$f")
            plan_indexes[$idx]=$(b64decode "$g")
            ;;
        "")
            ;;
        *)
            echo "ERROR: invalid BLUE bundle record: $kind"
            exit 1
            ;;
    esac
done < "$bundle"

[[ $population_meta == "$pop" ]] || {
    echo "ERROR: bundle population mismatch"
    exit 1
}

mapfile -t file_indices < <(
    printf '%s\n' "${!file_name[@]}" | sort -n
)
mapfile -t plan_indices < <(
    printf '%s\n' "${!plan_id[@]}" | sort -n
)

((${#file_indices[@]} >= 2)) || {
    echo "ERROR: BLUE requires >=2 files"
    exit 1
}
((${#plan_indices[@]} >= 1)) || {
    echo "ERROR: no selected BLUE plan"
    exit 1
}
((${#traits[@]} >= 1)) || {
    echo "ERROR: no cleaned traits"
    exit 1
}

# Selected source paths from web page.
[[ -s "$staging/data-source.txt" ]] || {
    echo "ERROR: missing data-source.txt"
    exit 1
}
mapfile -t source_fields < "$staging/data-source.txt"
((${#source_fields[@]} == 5)) || {
    echo "ERROR: invalid data-source.txt"
    exit 1
}
sample_source=${source_fields[1]}
[[ -s "$sample_source" ]] || {
    echo "ERROR: missing sampleList: $sample_source"
    exit 1
}

# Project program folder follows upstream behavior.
canonical_program="$BASE/000data_prepare/program"
[[ -d "$canonical_program" ]] || {
    echo "ERROR: missing canonical program: $canonical_program"
    exit 1
}
if [[ ! -e "$work/program" ]]; then
    cp -aL "$canonical_program" "$work/program"
    echo "PROJECT_PROGRAM_COPY: $canonical_program -> $work/program"
fi

# BLUE R scripts always come from the canonical server directory.
blue_program_dir="/data9/home/yzhao/program/EMMAx"
blue_multi_env="$blue_program_dir/generate_blue_multi_env_allTrait.R"
blue_locyear="$blue_program_dir/generate_blue_multiEnvYear_allTrait.R"

[[ -s "$blue_multi_env" ]] || {
    echo "ERROR: missing BLUE R program: $blue_multi_env"
    exit 1
}
[[ -s "$blue_locyear" ]] || {
    echo "ERROR: missing BLUE R program: $blue_locyear"
    exit 1
}

cp -f "$blue_multi_env" \
    "$work/program/generate_blue_multi_env_allTrait.R"
cp -f "$blue_locyear" \
    "$work/program/generate_blue_multiEnvYear_allTrait.R"

echo "BLUE_PROGRAM_SOURCE: $blue_program_dir"
echo "BLUE_PROGRAM_READY: $work/program/generate_blue_multi_env_allTrait.R"
echo "BLUE_PROGRAM_READY: $work/program/generate_blue_multiEnvYear_allTrait.R"

RSCRIPT=/public/home/yzhao/tool/R-4.1.2/bin/Rscript
[[ -x "$RSCRIPT" ]] || {
    echo "ERROR: missing Rscript: $RSCRIPT"
    exit 1
}

# Direct lme4 conversion.
# Only the environment prefix is removed from trait headers.
cat > "$staging/make_lme4_input.R" <<'RS'
args <- commandArgs(TRUE)
infile <- args[1]
outfile <- args[2]
prefix <- args[3]

first <- readLines(infile, n=1, warn=FALSE)
if (!length(first)) stop("Empty phenotype file")
sep <- if (grepl("	", first, fixed=TRUE)) "	" else ","

d <- read.table(
    infile,
    header=TRUE,
    sep=sep,
    check.names=FALSE,
    quote="\"",
    comment.char="",
    stringsAsFactors=FALSE,
    na.strings=c("NA","")
)

if (ncol(d) < 2 || names(d)[1] != "Accession") {
    stop("First column must be Accession")
}

nm <- names(d)

for (i in 2:length(nm)) {
    if (nzchar(prefix) && startsWith(nm[i], prefix)) {
        nm[i] <- substring(nm[i], nchar(prefix) + 1L)
    }
}

if (any(!nzchar(nm[-1]))) stop("Empty trait name after prefix removal")
if (anyDuplicated(nm[-1])) stop("Duplicated trait names after prefix removal")

names(d) <- nm
write.csv(d, outfile, row.names=FALSE, quote=FALSE, na="NA")
RS

declare -A lme4_file

for idx in "${file_indices[@]}"; do
    src="$staging/blue-unpack/$idx.raw"
    name=${file_name[$idx]}
    prefix=${file_prefix[$idx]}
    loc=${file_loc[$idx]}
    year=${file_year[$idx]}

    [[ -s "$src" ]] || {
        echo "ERROR: empty unpacked input: $name"
        exit 1
    }

    # Publish original uploaded phenotype.
    target="$work/phenotype/$name"
    [[ ! -d "$target" ]] || {
        echo "ERROR: refusing to replace directory: $target"
        exit 1
    }
    [[ ! -L "$target" ]] || rm -f "$target"
    cp -f "$src" "$target"
    cmp -s "$src" "$target" || {
        echo "ERROR: phenotype publish verification failed: $name"
        exit 1
    }
    echo "UPLOAD_READY: $target"

    stem=${name%.*}
    # Ensure lme4 filenames always contain current population + LocYear so the
    # LocYear R program never needs a hard-coded GPallInd prefix.
    if [[ $stem == "${pop}_"* ]]; then
        lme4_name="${stem}_lme4.csv"
    else
        lme4_name="${pop}_${stem}_lme4.csv"
    fi

    out="$work/lme4/$lme4_name"
    echo "LME4_CONVERT: $name -> $lme4_name"
    echo "  prefix=$prefix"
    echo "  environment=${loc}${year}"

    "$RSCRIPT" "$staging/make_lme4_input.R" \
        "$target" "$out" "$prefix"

    [[ -s "$out" ]] || {
        echo "ERROR: empty lme4 output: $out"
        exit 1
    }

    lme4_file[$idx]=$out
done

# Verify every lme4 file has the identical cleaned header.
reference_header=""
for idx in "${file_indices[@]}"; do
    h=$(head -n 1 "${lme4_file[$idx]}" | tr -d '\r')
    if [[ -z $reference_header ]]; then
        reference_header=$h
    elif [[ $h != "$reference_header" ]]; then
        echo "ERROR: cleaned lme4 headers differ"
        echo "REFERENCE: $reference_header"
        echo "CURRENT:   $h"
        exit 1
    fi
done

# Helper: normalize one BLUE raw result for combined GWAS.
# Historical R output is:
#   <trait>_<plan-label>
# Combined GWAS needs:
#   <plan-label>_<trait>
cat > "$staging/normalize_blue_for_gwas.R" <<'RS'
args <- commandArgs(TRUE)
infile <- args[1]
outfile <- args[2]
label <- args[3]

d <- read.delim(
    infile,
    header=TRUE,
    check.names=FALSE,
    stringsAsFactors=FALSE,
    quote="",
    comment.char=""
)

if (ncol(d) < 2) stop("BLUE file contains no trait columns")
suffix <- paste0("_", label)
nm <- names(d)

for (i in 2:length(nm)) {
    trait <- nm[i]
    if (endsWith(trait, suffix)) {
        trait <- substr(trait, 1L, nchar(trait) - nchar(suffix))
    }
    nm[i] <- paste(label, trait, sep="_")
}

if (anyDuplicated(nm[-1])) stop("Duplicated normalized BLUE trait names")
names(d) <- nm

write.table(
    d,
    outfile,
    row.names=FALSE,
    quote=FALSE,
    sep="\t",
    na="NA"
)
RS

normalized_files=()

for pidx in "${plan_indices[@]}"; do
    kind=${plan_kind[$pidx]}
    label=${plan_label[$pidx]}
    output=${plan_output[$pidx]}
    title=${plan_title[$pidx]}
    index_csv=${plan_indexes[$pidx]}

    IFS=',' read -r -a indexes <<< "$index_csv"
    input_files=()
    for idx in "${indexes[@]}"; do
        [[ -n ${lme4_file[$idx]+x} ]] || {
            echo "ERROR: plan references missing file index $idx"
            exit 1
        }
        input_files+=("${lme4_file[$idx]}")
    done

    raw="$work/lme4/blue/${output}.txt"
    filtered="$work/lme4/blue/${output}.csv"
    normalized="$staging/${output}.gwas.txt"

    echo "============================================================"
    echo "BLUE_PLAN: $title"
    echo "BLUE_KIND: $kind"
    echo "BLUE_LABEL: $label"
    echo "BLUE_OUTPUT: $raw"
    printf 'BLUE_INPUT: %s\n' "${input_files[@]}"

    if [[ $kind == multi_env ]]; then
        echo "BLUE_MODEL: trait ~ Accession + (1 | env)"
        "$RSCRIPT" \
            "$work/program/generate_blue_multi_env_allTrait.R" \
            "$raw" "$label" "${input_files[@]}"

    elif [[ $kind == locyear ]]; then
        echo "BLUE_MODEL: trait ~ Accession + (1 | Loc) + (1 | Year) + (1 | Loc:Year)"
        "$RSCRIPT" \
            "$work/program/generate_blue_multiEnvYear_allTrait.R" \
            "$raw" "$label" "${input_files[@]}"

    else
        echo "ERROR: unknown BLUE plan kind: $kind"
        exit 1
    fi

    [[ -s "$raw" ]] || {
        echo "ERROR: BLUE output is empty: $raw"
        exit 1
    }
    echo "BLUE_TXT_READY: $raw"

    # Historical sample-filtered CSV retained exactly as before.
    perl "$work/program/select.pheno.only.pl" \
        "$sample_source" "$raw" "$filtered"
    [[ -s "$filtered" ]] || {
        echo "ERROR: BLUE sample-filtered CSV is empty: $filtered"
        exit 1
    }
    echo "BLUE_CSV_READY: $filtered"

    "$RSCRIPT" "$staging/normalize_blue_for_gwas.R" \
        "$raw" "$normalized" "$label"
    [[ -s "$normalized" ]] || {
        echo "ERROR: normalized BLUE GWAS file is empty"
        exit 1
    }
    normalized_files+=("$normalized")
done

# Merge all selected plans by Accession.
cat > "$staging/merge_blue_for_gwas.R" <<'RS'
args <- commandArgs(TRUE)
outfile <- args[1]
files <- args[-1]

if (length(files) < 1) stop("No normalized BLUE files")

all <- NULL
for (f in files) {
    d <- read.delim(
        f,
        header=TRUE,
        check.names=FALSE,
        stringsAsFactors=FALSE,
        quote="",
        comment.char=""
    )
    if (!("Accession" %in% names(d))) {
        stop(paste("Missing Accession in", f))
    }

    if (is.null(all)) {
        all <- d
    } else {
        all <- merge(all, d, by="Accession", all=TRUE, sort=FALSE)
    }
}

if (anyDuplicated(names(all)[-1])) stop("Duplicate traits after BLUE merge")
all <- all[order(all$Accession), , drop=FALSE]

write.table(
    all,
    outfile,
    row.names=FALSE,
    quote=FALSE,
    sep="\t",
    na="NA"
)
RS

combined="$staging/${plan_output[${plan_indices[0]}]}.txt"
"$RSCRIPT" "$staging/merge_blue_for_gwas.R" \
    "$combined" "${normalized_files[@]}"

[[ -s "$combined" ]] || {
    echo "ERROR: combined BLUE phenotype is empty"
    exit 1
}

echo "============================================================"
echo "BLUE_COMBINED_GWAS_READY: $combined"
echo "BLUE_SELECTED_PLANS: ${#plan_indices[@]}"
echo "BLUE_TOTAL_GWAS_TRAITS: $(( ${#plan_indices[@]} * ${#traits[@]} ))"
echo "============================================================"

# Preserve original bundle, then replace staging phenotype with combined BLUE.
cp -f "$bundle" "$staging/blue-bundle.txt"
cp -f "$combined" "$staging/phenotype.upload"
cmp -s "$combined" "$staging/phenotype.upload" || {
    echo "ERROR: failed to stage combined BLUE phenotype"
    exit 1
}

# Audit configuration.
{
    echo "mode=multi-blue"
    echo "population=$pop"
    echo "inputCount=${#file_indices[@]}"
    echo "planCount=${#plan_indices[@]}"
    for idx in "${file_indices[@]}"; do
        printf 'input=%s	loc=%s	year=%s	prefix=%s
' \
            "${file_name[$idx]}" \
            "${file_loc[$idx]}" \
            "${file_year[$idx]}" \
            "${file_prefix[$idx]}"
    done
    for pidx in "${plan_indices[@]}"; do
        printf 'plan=%s\tkind=%s\tlabel=%s\toutput=%s\n' \
            "${plan_id[$pidx]}" \
            "${plan_kind[$pidx]}" \
            "${plan_label[$pidx]}" \
            "${plan_output[$pidx]}"
    done
} > "$staging/blue-config.txt"

echo "BLUE_PREP_COMPLETE"
echo "DELEGATE_TO_ORIGINAL_PREPARE"

# Original workflow now owns prepare.exit and downstream state.
trap - EXIT
exec bash "$ORIGINAL" "$@"
