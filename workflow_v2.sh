#!/usr/bin/env bash
set -Eeuo pipefail
export LC_ALL=C
BASE=/data9/home/yzhao/GWAS_IRGSP1.0
stage=$1; pop=$2; folder=$3; filename=$4; runid=$5
[[ $pop =~ ^(GPall|GPallInd|GPallJap)$ && $folder =~ ^[A-Za-z0-9][A-Za-z0-9_.-]{0,79}$ && $runid =~ ^[A-Za-z0-9_-]+$ ]] || exit 64
stem=${filename%.*}
[[ $stem =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || exit 64
tag=${6:-$stem}
[[ $tag =~ ^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$ ]] || exit 64
work="$BASE/$folder"
staging="$BASE/.gwas-web/$runid"
batch="$work/.gwas-runs/$runid"
manifest="$batch/manifest.tsv"
archive="${stem}__${runid}"
if [[ -s $staging/archive-name.txt ]]; then
  archive=$(cat "$staging/archive-name.txt")
  [[ $archive =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*_[0-9]{8}-[0-9]{4}$ ]] || exit 64
fi
lock="$work/.gwas-active"
mkdir -p "$staging"
exec >> "$staging/run.log" 2>&1
PS4='+ ${LINENO}: '
echo "============================================================"
echo "STAGE_START: $stage"
echo "RUN_ID: $runid"
echo "POPULATION: $pop"
echo "PHENOTYPE: $filename"
echo "PROJECT_DIR: $work"
echo "============================================================"
die() { echo "ERROR: $*"; exit 1; }
finish() {
  rc=$?
  trap - EXIT
  printf '%s\n' "$rc" > "$staging/$stage.exit.tmp"
  mv "$staging/$stage.exit.tmp" "$staging/$stage.exit"
  if ((rc != 0)); then echo "STOP: run=$runid; inspect $batch and queue before unlocking $lock"; fi
}
trap finish EXIT
command -v perl >/dev/null
[[ $(id -un) == yzhao ]] || die 'Expected remote user yzhao'
shopt -s nullglob
src="$BASE/000data_prepare/EMMAx.Data/${pop}_miss20"
source_version="${pop}miss20-864lines"
sample_source="$src/$pop.sampleList"
tped_template="$src/${pop}_chr{chr}.tped"
kin_source="$src/${pop}_allChr_snpNonHet.hIBS.kinf"
maf_template="/public/home/yzhao/IRGSP1.0_GPall/MAF/${pop}_chr{chr}_MAF.vcf"
if [[ -f $staging/data-source.txt ]]; then
  mapfile -t source_fields < "$staging/data-source.txt"
  ((${#source_fields[@]} == 5)) || die 'Invalid data source configuration'
  source_version=${source_fields[0]}; sample_source=${source_fields[1]}
  tped_template=${source_fields[2]}; kin_source=${source_fields[3]}; maf_template=${source_fields[4]}
fi
source_path() { local template=$1; printf '%s' "${template//\{chr\}/$2}"; }
if [[ $stage == inflation ]]; then
  work="$BASE/$folder"
  archive_dir="$work/EMMAx.Result/hIBS/$stem"
  output_dir="$archive_dir/output"
  traits_file="$archive_dir/traits.txt"
  canonical_inflation_script="$BASE/000data_prepare/program/get_genomicInflationFactor.py"
  inflation_program_dir="$work/.gwas-runs/$runid/program"
  inflation_script="$inflation_program_dir/get_genomicInflationFactor.py"
  inflation_python="/data6/tool/anaconda2-4.1.1/bin/python"
  inflation_dir="$work/inflationFactor"

  echo "============================================================"
  echo "INFLATION_ONLY_START"
  echo "PROJECT_DIR: $work"
  echo "RESULT_ARCHIVE: $archive_dir"
  echo "INFLATION_DIR: $inflation_dir"
  echo "INFLATION_PYTHON: $inflation_python"
  echo "============================================================"

  [[ -d "$output_dir" ]] || die "Missing archived EMMAX output directory: $output_dir"
  [[ -s "$traits_file" ]] || die "Missing archived trait list: $traits_file"
  [[ -s "$canonical_inflation_script" ]] \
    || die "Missing genomic inflation program: $canonical_inflation_script"
  [[ -x "$inflation_python" ]] \
    || die "Missing inflation Python: $inflation_python"

  mkdir -p "$inflation_dir" "$inflation_program_dir"
  cp -f "$canonical_inflation_script" "$inflation_script"
  [[ -s "$inflation_script" ]] \
    || die "Failed to copy genomic inflation program to: $inflation_script"

  echo "INFLATION_PROGRAM_SOURCE: $canonical_inflation_script"
  echo "INFLATION_PROGRAM_READY: $inflation_script"

  "$inflation_python" -c 'import numpy, scipy' \
    || die "Inflation Python is missing numpy/scipy: $inflation_python"

  while read -r trait; do
    ps_pattern="$output_dir/${trait}_chr01_hIBS.ps"
    inflation_out="$inflation_dir/${trait}.inflationFactor"

    echo "INFLATION_FACTOR_START: $trait"
    echo "INFLATION_FACTOR_INPUT: $ps_pattern"
    echo "INFLATION_FACTOR_OUTPUT: $inflation_out"

    [[ -s "$ps_pattern" ]] \
      || die "Missing chr01 EMMAX ps file for $trait: $ps_pattern"

    "$inflation_python" "$inflation_script" "$ps_pattern" "$inflation_out" \
      || die "Genomic inflation calculation failed for trait: $trait"

    [[ -s "$inflation_out" ]] \
      || die "Empty inflation factor output for trait: $trait"
    echo "INFLATION_FACTOR_DONE: $trait"
  done < "$traits_file"

  echo "INFLATION_FACTOR_ALL_DONE: $inflation_dir"
  echo "SUCCESS: inflation factor calculation completed from archived EMMAX ps files"
  exit 0
fi

if [[ $stage == prepare ]]; then
  mkdir -p "$work"
  mkdir "$lock" || die 'Another run or unresolved failure owns this directory; inspect .gwas-active'
  printf '%s\n' "$runid" > "$lock/owner"
  [[ ! -e $batch ]] || die 'Run ID already exists; will not replay'
  mkdir -p "$batch"
  cd "$work"
  mkdir -p chromosome EMMAx.Data EMMAx.Result/hIBS/output phenotype plot QQplot sh.file sigSNP tped
  if [[ -f .gwas-population ]]; then
    [[ $(cat .gwas-population) == "$pop" ]] || die 'Existing directory uses another population'
  fi
  # Publish and verify the phenotype before any large TPED/kinship copies.
  # The web layer has already handled the one-time same-name confirmation.
  if [[ -L phenotype/$filename ]]; then rm -f "phenotype/$filename"; fi
  [[ ! -d phenotype/$filename ]] || die "Refusing to replace a directory: phenotype/$filename"
  cp -f "$staging/phenotype.upload" "phenotype/$filename"
  [[ -s phenotype/$filename ]] || die 'Uploaded phenotype is empty'
  cmp -s "$staging/phenotype.upload" "phenotype/$filename" || die 'Uploaded phenotype verification failed'
  echo "UPLOAD_READY: $work/phenotype/$filename"
  printf '%s\n' "$source_version" "$sample_source" "$tped_template" "$kin_source" "$maf_template" > "$batch/data-source.txt"
  sources=("$sample_source" "$kin_source")
  for chr in {01..12}; do sources+=("$(source_path "$tped_template" "$chr")" "$(source_path "$maf_template" "$chr")"); done
  for f in "${sources[@]}"; do [[ -s $f ]] || die "Missing or empty data source: $f"; done
  # User-owned analysis programs always come from the canonical server directory.
  # Nothing under program/ is uploaded from the Windows web console.
  canonical_program="$BASE/000data_prepare/program"
  [[ -d $canonical_program ]] || die "Missing canonical server program directory: $canonical_program"
  if [[ ! -e program ]]; then
    cp -aL "$canonical_program" program
    echo "PROJECT_PROGRAM_COPY: $canonical_program -> $work/program"
  fi
  cp -aL "$canonical_program" "$batch/program"
  echo "BATCH_PROGRAM_SOURCE: $canonical_program"
  checked_copy() {
    local from=$1 to=$2
    [[ ! -d $to ]] || die "Input destination is a directory: $to"
    cp "$from" "$to.tmp.$runid"
    cmp -s "$from" "$to.tmp.$runid" || die "Copy verification failed: $to"
    mv -T "$to.tmp.$runid" "$to"
  }
  checked_copy "$sample_source" "phenotype/$pop.sampleList"
  checked_copy "$kin_source" "tped/${pop}_allChr_snpNonHet.hIBS.kinf"
  for chr in {01..12}; do checked_copy "$(source_path "$tped_template" "$chr")" "tped/${pop}_chr${chr}.tped"; done
  echo "INPUTS_READY:"
  echo "  sampleList: $sample_source"
  echo "  kinship:    $kin_source"
  echo "  TPED:       ${tped_template//\{chr\}/01} ... ${tped_template//\{chr\}/12} (12 chromosomes)"
  cp "$batch/data-source.txt" .gwas-data-source.txt
  printf '%s\n' "$pop" > .gwas-population
  for tool in select.pheno.only.pl EMMAx-Step2.ChangetoEMMAxPhenotype.tfam.pl qsub_nodeAdmin.pl EMMAx.run-get.pValue.hIBS.pl plot_Allpicture_ofOneTrait.pl plot_manhattan.pl plot_manhattan_eachChr.pl plot_QQplot2023.pl stat.allSNPSignificantPos.inOneDir.list.pl get_genomicInflationFactor.py; do
    [[ -s $batch/program/$tool ]] || die "Missing batch program/$tool"
  done
  # Empty conversion workspace isolates this upload, including duplicate trait headers.
  convert="$batch/convert"
  mkdir -p "$convert"/{EMMAx.Data,EMMAx.Result,phenotype,chromosome,plot,QQplot,sh.file,sigSNP}
  ln -s "$batch/program" "$convert/program"
  ln -s "$work/tped" "$convert/tped"
  ln -s "$work/phenotype/$pop.sampleList" "$convert/phenotype/$pop.sampleList"
  ln -s "$work/phenotype/$filename" "$convert/phenotype/$filename"
  cd "$convert"
  echo "COMMAND: perl program/select.pheno.only.pl phenotype/$pop.sampleList phenotype/$filename phenotype/${pop}_${stem}.csv"
  perl program/select.pheno.only.pl "phenotype/$pop.sampleList" "phenotype/$filename" "phenotype/${pop}_${stem}.csv"
  [[ -s phenotype/${pop}_${stem}.csv ]] || die 'Empty extracted phenotype'
  echo "COMMAND: perl program/EMMAx-Step2.ChangetoEMMAxPhenotype.tfam.pl phenotype/${pop}_${stem}.csv $pop"
  perl program/EMMAx-Step2.ChangetoEMMAxPhenotype.tfam.pl "phenotype/${pop}_${stem}.csv" "$pop"
  generated=(EMMAx.Data/"${pop}_"*.tped)
  ((${#generated[@]})) || die 'No TPED generated for this upload'
  : > "$batch/raw.tsv"
  for f in "${generated[@]}"; do
    name=${f##*/}; name=${name%.tped}
    [[ $name =~ ^(${pop}_[A-Za-z0-9_.-]+)_chr([0-9]{1,2})$ ]] || die "Unexpected TPED name $name"
    old=${BASH_REMATCH[1]}; chr=${BASH_REMATCH[2]}
    chr=$((10#$chr)); ((chr >= 1 && chr <= 12)) || die "Invalid rice chromosome: $chr"
    printf -v chr '%02d' "$chr"
    if [[ $old == "${pop}_${tag}" || $old == "${pop}_${tag}_"* ]]; then trait=$old
    else trait="${pop}_${tag}_${old#${pop}_}"; fi
    printf '%s\t%s\t%s\t%s\n' "$old" "$trait" "$chr" "$name" >> "$batch/raw.tsv"
  done
  sort -t $'\t' -k2,2 -k3,3 "$batch/raw.tsv" > "$manifest"
  cut -f2 "$manifest" | sort -u > "$batch/traits.txt"
  [[ -s $staging/expected_traits.txt ]] || die 'Missing expected upload trait list'
  sort "$staging/expected_traits.txt" > "$batch/expected.sorted"
  diff -u "$batch/expected.sorted" "$batch/traits.txt" || die 'Converted traits do not exactly match the uploaded phenotype columns'
  while read -r trait; do
    count=$(awk -F '\t' -v t="$trait" '$2==t {n++} END {print n+0}' "$manifest")
    [[ $count == 12 ]] || die "Trait $trait requires exactly 12 distinct chromosomes; found $count"
    for chr in {01..12}; do
      count=$(awk -F '\t' -v t="$trait" -v c="$chr" '$2==t && $3==c {n++} END {print n+0}' "$manifest")
      [[ $count == 1 ]] || die "Trait $trait missing or duplicated chr$chr"
    done
  done < "$batch/traits.txt"
  # Validate generated inputs, but do not stop merely because an exact target
  # from an earlier run exists.  The user has already confirmed same-name rerun.
  # No directory-wide cleanup is performed.
  while IFS=$'\t' read -r old trait chr raw; do
    [[ -s EMMAx.Data/$raw.tped ]] || die "Empty TPED $raw"
    [[ -s EMMAx.Data/$raw.tfam && -s EMMAx.Data/$old.tfam ]] || die "Missing TFAM for $raw"
    [[ -s EMMAx.Data/$old.pheno || -s $old.pheno ]] || die "Missing phenotype $old.pheno"
  done < "$manifest"
  while IFS=$'\t' read -r old trait chr raw; do
    for f in EMMAx.Data/"$raw".*; do
      ext=${f##*.}
      target="$work/EMMAx.Data/${trait}_chr${chr}.$ext"
      [[ ! -d $target || -L $target ]] || die "Refusing to replace directory target: $target"
      rm -f "$target"
      ln -s "$convert/$f" "$target"
    done
  done < "$manifest"
  while read -r trait; do
    old=$(awk -F '\t' -v t="$trait" '$2==t {print $1; exit}' "$manifest")
    f="EMMAx.Data/$old.pheno"; [[ -s $f ]] || f="$old.pheno"
    for target in "$work/EMMAx.Data/$trait.pheno" "$work/EMMAx.Data/$trait.tfam"; do
      [[ ! -d $target || -L $target ]] || die "Refusing to replace directory target: $target"
      rm -f "$target"
    done
    ln -s "$convert/$f" "$work/EMMAx.Data/$trait.pheno"
    ln -s "$convert/EMMAx.Data/$old.tfam" "$work/EMMAx.Data/$trait.tfam"
  done < "$batch/traits.txt"
  cp -f "phenotype/${pop}_${stem}.csv" "$work/phenotype/"
  echo "MANIFEST_SAVED: $manifest"
  echo "TRAITS:"
  sed 's/^/  /' "$batch/traits.txt"
  printf 'BATCH_READY: phenotype=%s; population=%s; traits=%s; chromosome_jobs=%s\n' "$filename" "$pop" "$(wc -l < "$batch/traits.txt")" "$(wc -l < "$manifest")"
  exit 0
fi
[[ -f $lock/owner && $(cat "$lock/owner") == "$runid" ]] || die 'Run does not own the directory lock'
[[ -s $manifest && -s $batch/traits.txt ]] || die 'Missing frozen manifest; refusing a directory-wide scan'
cd "$work"
if [[ $stage == admin ]]; then
  [[ $(cat "$staging/prepare.exit") == 0 ]] || die 'Preparation incomplete'
  command -v qstat >/dev/null
  command -v qsub >/dev/null
  [[ ! -e $batch/admin.started ]] || die 'Admin phase already started; refusing duplicate submission'
  : > "$batch/admin.started"
  [[ ! -d sh.file/emmax/$archive && ! -d sh.file/plotGWAS/$archive ]] || die 'Queue-log archive name already exists; choose a new run time'
  # Existing phenotype-specific result/sigSNP targets are allowed after the
  # one-time web confirmation. They are replaced only when this run writes them.
  mkdir -p "$batch/status"
  cat > "$batch/worker.sh" <<'WORKER'
#!/usr/bin/env bash
set -Eeuo pipefail
cd "$1"
status=$2; key=$3
shift 3
done_marker() {
  rc=$?
  trap - EXIT
  printf '%s\n' "$rc" > "$status/$key.exit.tmp"
  mv "$status/$key.exit.tmp" "$status/$key.exit"
}
trap done_marker EXIT
"$@"
WORKER
  submit() {
    local key=$1 name=$2 command output id
    shift 2
    printf -v command '%q ' bash "$batch/worker.sh" "$work" "$batch/status" "$key" "$@"
    printf 'RUN_COMMAND:'
    printf ' %q' "$@"
    printf '\n'
    output=$(perl "$batch/program/qsub_nodeAdmin.pl" "$name" "$command") || die "Submission failed: $name"
    printf '%s\n' "$output" | tee -a "$batch/submissions.log"
    id=$(printf '%s\n' "$output" | sed -nE 's/^([0-9]+)(\.[A-Za-z0-9_.-]+)?[[:space:]]*$/\1/p; s/.*Your job(-array)? ([0-9]+).*/\2/p' | sort -u)
    [[ $id =~ ^[0-9]+$ ]] || die "Cannot identify exactly one queue job ID for $name; inspect queue before retrying: $output"
    printf '%s\t%s\t%s\n' "$key" "$id" "$name" >> "$batch/$phase.jobs"
    echo "TRACK_JOB: phase=$phase id=$id script=$name"
    sleep 1
  }
  phase_logs() {
    logs=()
    while read -r key id name; do
      found=("$name".o* "$name".e*)
      ((${#found[@]})) || die "No scheduler logs for $name"
      logs+=("${found[@]}")
    done < "$batch/$phase.jobs"
  }
  wait_jobs() {
    local start=$SECONDS snapshot key id name code state active incomplete completed total missing_since=0 wait_seconds
    total=$(wc -l < "$batch/$phase.jobs")
    while true; do
      snapshot=$(qstat -u yzhao) || die 'qstat failed; refusing to infer completion'
      printf '%s\n' "$snapshot" > "$batch/queue.latest"

      active=0
      incomplete=0
      completed=0

      while read -r key id name; do
        state=$(printf '%s\n' "$snapshot" | awk -v id="$id" '
          {
            split($1,a,".")
            if (a[1]==id) {
              for (i=4;i<=NF;i++) {
                if ($i ~ /^[A-Z]$/) s=$i
              }
              print s
              exit
            }
          }')

        # qstat "C" is useful information, but is NOT required for completion.
        # Some clusters remove completed jobs before the next 120 s poll.
        if [[ -n $state && $state != C ]]; then
          active=$((active + 1))
        fi

        if [[ -f $batch/status/$key.exit ]]; then
          code=$(tr -d '\r\n[:space:]' < "$batch/status/$key.exit")
          [[ $code =~ ^-?[0-9]+$ ]] || die "Invalid completion marker for job $id ($key): $code"
          if [[ $code != 0 ]]; then
            found=("$name".o* "$name".e*)
            if ((${#found[@]})); then tail -n 100 "${found[@]}"; fi
            die "Job $id ($key) failed with exit code $code"
          fi
          completed=$((completed + 1))
        else
          incomplete=$((incomplete + 1))
        fi
      done < "$batch/$phase.jobs"

      # Primary success rule:
      # every worker wrote exit=0 AND this batch has no non-C job still in qstat.
      if ((incomplete == 0 && active == 0)); then
        echo "QUEUE_PHASE_COMPLETE: $phase; completed=$completed/$total"
        break
      fi

      # If jobs have vanished from qstat, allow a short race window for worker.sh
      # to atomically write its .exit marker. Do not wait for qstat C.
      if ((active == 0 && incomplete > 0)); then
        if ((missing_since == 0)); then
          missing_since=$SECONDS
        fi
        if ((SECONDS - missing_since >= 60)); then
          die "Jobs left qstat but $incomplete/$total completion markers are still missing after 60 seconds"
        fi
        wait_seconds=10
      else
        missing_since=0
        wait_seconds=120
      fi

      ((SECONDS-start < 604800)) || die 'Monitoring exceeded 7 days; jobs are not cancelled'
      echo "WAIT: $phase; completed=$completed/$total; active_in_qstat=$active; missing_markers=$incomplete; next_check=${wait_seconds}s"
      sleep "$wait_seconds"
    done

    # Give scheduler stdout/stderr a moment to finish flushing.
    sleep 1
    phase_logs
    if grep -inE 'error|fatal|segmentation fault|traceback' "${logs[@]}"; then
      die "Error text detected after $phase; this batch logs are preserved"
    else
      code=$?; [[ $code == 1 ]] || die 'Failed to read scheduler logs'
    fi
  }
  phase=emmax
  echo '===== GWAS SUBMISSION / MONITORING ====='
  while IFS=$'\t' read -r old trait chr raw; do
    [[ -s EMMAx.Data/${trait}_chr${chr}.tped && -s EMMAx.Data/$trait.pheno ]] || die "Missing batch input $trait chr$chr"
    [[ -s $(source_path "$maf_template" "$chr") ]] || die "Missing MAF chr$chr"
  done < "$manifest"
  : > "$batch/emmax.jobs"
  while IFS=$'\t' read -r old trait chr raw; do
    name="${trait}_chr${chr}"
    submit "${name}_emmax" "${name}_emmax_hibs.sh" perl "$batch/program/EMMAx.run-get.pValue.hIBS.pl" "$name" "$trait.pheno" "tped/${pop}_allChr_snpNonHet.hIBS.kinf" "$(source_path "$maf_template" "$chr")"
  done < "$manifest"
  wait_jobs
  while IFS=$'\t' read -r old trait chr raw; do
    [[ -s EMMAx.Result/hIBS/${trait}_chr${chr}.hIBS ]] || die "Missing result $trait chr$chr"
  done < "$manifest"
  mkdir -p "sh.file/emmax/$archive"
  mv "${logs[@]}" "sh.file/emmax/$archive/"
  while read -r key id name; do [[ ! -f $name ]] || mv "$name" "sh.file/emmax/$archive/"; done < "$batch/emmax.jobs"
  phase=plot
  echo '===== PLOT SUBMISSION / MONITORING ====='
  # The supplied umbrella script ignores system() failures. Propagate each one.
  cat > "$batch/plot_checked.sh" <<'PLOT'
#!/usr/bin/env bash
set -Eeuo pipefail
program_dir=$1
shift
for script in plot_manhattan.pl plot_manhattan_eachChr.pl plot_QQplot2023.pl; do
  perl "$program_dir/$script" "$@"
done
PLOT
  : > "$batch/plot.jobs"
  while read -r trait; do
    submit "${trait}_plot" "${trait}_plotGWAS.sh" bash "$batch/plot_checked.sh" "$batch/program" "EMMAx.Result/hIBS/${trait}_chr01.hIBS" "$trait" "$trait"
  done < "$batch/traits.txt"
  wait_jobs
  while read -r trait; do
    [[ -s chromosome/$trait.png && -s QQplot/${trait}_QQplot.png && -s QQplot/$trait.value ]] || die "Missing Manhattan/QQ outputs for $trait"
    for chr in {01..12}; do
      [[ -s plot/${trait}_chr${chr}.png ]] || die "Missing chromosome plot $trait chr$chr"
    done
  done < "$batch/traits.txt"
  mkdir -p "sh.file/plotGWAS/$archive"
  mv "${logs[@]}" "sh.file/plotGWAS/$archive/"
  while read -r key id name; do [[ ! -f $name ]] || mv "$name" "sh.file/plotGWAS/$archive/"; done < "$batch/plot.jobs"
elif [[ $stage == final ]]; then
  [[ $(cat "$staging/admin.exit") == 0 ]] || die 'Queue/plot phase incomplete'
  [[ ! -e $batch/final.started ]] || die 'Final phase already started; inspect outputs before recovery'
  : > "$batch/final.started"
  final="$batch/final"
  mkdir -p "$final/EMMAx.Result/hIBS" "$final/sigSNP"
  ln -s "$batch/program" "$final/program"
  while IFS=$'\t' read -r old trait chr raw; do
    name="${trait}_chr${chr}.hIBS"
    ln -s "$work/EMMAx.Result/hIBS/$name" "$final/EMMAx.Result/hIBS/$name"
  done < "$manifest"
  cd "$final"
  peak_prefix="${pop}_${stem}"

  echo "COMMAND: perl program/stat.allSNPSignificantPos.inOneDir.list.pl EMMAx.Result/hIBS hIBS $pop $pop"
  perl program/stat.allSNPSignificantPos.inOneDir.list.pl EMMAx.Result/hIBS hIBS "$pop" "$pop"
  [[ -f sigSNP/$pop.list ]] || die 'Missing peak SNP list'

  # The original Perl program writes <population>.list.
  # Keep the original program unchanged, but rename this batch's output to
  # <population>_<phenotype>.list for traceable multi-phenotype archives.
  mv "sigSNP/$pop.list" "sigSNP/${peak_prefix}.list"
  echo "PEAK_LIST: sigSNP/${peak_prefix}.list"

  echo "COMMAND: perl /data5/home/yzhao/program/compared_peakSNP_with_RiceNaviGene_1.0.pl sigSNP/${peak_prefix}.list sigSNP/${peak_prefix}_geneList.csv"
  perl /data5/home/yzhao/program/compared_peakSNP_with_RiceNaviGene_1.0.pl     "sigSNP/${peak_prefix}.list"     "sigSNP/${peak_prefix}_geneList.csv"
  [[ -f sigSNP/${peak_prefix}_geneList.csv ]] || die 'Missing gene annotation output'
  echo "GENE_LIST: sigSNP/${peak_prefix}_geneList.csv"
  cd "$work"

  # Peak SNP outputs stay directly under sigSNP/. Replace only these exact
  # phenotype-specific files; never clean the whole sigSNP directory.
  [[ ! -d "sigSNP/${peak_prefix}.list" ]] || die "Refusing to replace directory: sigSNP/${peak_prefix}.list"
  [[ ! -d "sigSNP/${peak_prefix}_geneList.csv" ]] || die "Refusing to replace directory: sigSNP/${peak_prefix}_geneList.csv"
  cp -f "$final/sigSNP/${peak_prefix}.list" "sigSNP/${peak_prefix}.list"
  cp -f "$final/sigSNP/${peak_prefix}_geneList.csv" "sigSNP/${peak_prefix}_geneList.csv"

  # GWAS result archive uses only the phenotype filename stem, e.g. LSH2022.
  dest="$work/EMMAx.Result/hIBS/$stem"
  [[ ! -e $dest || -d $dest ]] || die "Refusing to replace non-directory result archive: $dest"
  mkdir -p "$dest/output"
  while IFS=$'\t' read -r old trait chr raw; do
    name="${trait}_chr${chr}"
    mv "EMMAx.Result/hIBS/$name.hIBS" "$dest/"
    outputs=("EMMAx.Result/hIBS/output/${name}_hIBS."*)
    if ((${#outputs[@]})); then mv "${outputs[@]}" "$dest/output/"; fi
    ln -sfn "$dest/$name.hIBS" "$final/EMMAx.Result/hIBS/$name.hIBS"
  done < "$manifest"
  cp "$manifest" "$batch/traits.txt" "$dest/"

  echo '===== GENOMIC INFLATION FACTOR ====='
  inflation_script="$batch/program/get_genomicInflationFactor.py"
  [[ -s "$inflation_script" ]] || die "Missing genomic inflation program: $inflation_script"

  inflation_dir="$work/inflationFactor"
  mkdir -p "$inflation_dir"

  inflation_python="/data6/tool/anaconda2-4.1.1/bin/python"
  [[ -x "$inflation_python" ]] \
    || die "Missing inflation Python: $inflation_python"

  echo "INFLATION_PYTHON: $inflation_python"
  "$inflation_python" -c 'import numpy, scipy' \
    || die "Inflation Python is missing numpy/scipy: $inflation_python"

  while read -r trait; do
    ps_pattern="$dest/output/${trait}_chr01_hIBS.ps"
    inflation_out="$inflation_dir/${trait}.inflationFactor"

    echo "INFLATION_FACTOR_START: $trait"
    echo "INFLATION_FACTOR_INPUT: $ps_pattern"
    echo "INFLATION_FACTOR_OUTPUT: $inflation_out"

    [[ -s "$ps_pattern" ]] \
      || die "Missing chr01 EMMAX ps file for $trait: $ps_pattern"

    "$inflation_python" "$inflation_script" "$ps_pattern" "$inflation_out" \
      || die "Genomic inflation calculation failed for trait: $trait"

    [[ -s "$inflation_out" ]] \
      || die "Empty inflation factor output for trait: $trait"

    echo "INFLATION_FACTOR_DONE: $trait"
  done < "$batch/traits.txt"

  echo "INFLATION_FACTOR_ALL_DONE: $inflation_dir"

  # A recovered older batch may have archived queue logs under the old name.
  old_archive="${stem}__${runid}"
  for phase in emmax plot; do
    target_phase=$phase; [[ $phase != plot ]] || target_phase=plotGWAS
    target="sh.file/$target_phase/$archive"
    mkdir -p "$target"
    old_dir="sh.file/$phase/$old_archive"
    if [[ $old_dir != "$target" && -d $old_dir ]]; then
      old_files=("$old_dir/"*)
      for f in "${old_files[@]}"; do [[ ! -e $target/${f##*/} ]] || die "Archive collision: $f"; mv "$f" "$target/"; done
      rmdir "$old_dir"
    fi
    while read -r key id name; do
      if [[ -f $name ]]; then [[ ! -e $target/$name ]] || die "Script archive collision: $name"; mv "$name" "$target/"; fi
    done < "$batch/$phase.jobs"
  done
  printf 'RESULT_ARCHIVE: %s\nPEAK_LIST: %s\nGENE_LIST: %s\n' \
    "$dest" \
    "$work/sigSNP/${peak_prefix}.list" \
    "$work/sigSNP/${peak_prefix}_geneList.csv"
  # Only after all outputs pass validation, remove immediate input files.
  for input_dir in EMMAx.Data tped; do
    [[ ! -L $work/$input_dir && $(realpath "$work/$input_dir") == "$work/$input_dir" ]] || die "Unsafe cleanup directory: $input_dir"
    find "$work/$input_dir" -mindepth 1 -maxdepth 1 \( -type f -o -type l \) -print -delete
  done
  echo 'INPUT_CLEANUP_COMPLETE: EMMAx.Data and tped'
  mv "$lock/owner" "$batch/lock.completed"
  rmdir "$lock"
  echo 'SUCCESS: GWAS analysis completed; only this batch was processed and archived'
else die 'Unknown stage'; fi
