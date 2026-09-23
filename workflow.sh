#!/usr/bin/env bash
set -Eeuo pipefail
export LC_ALL=C
BASE=/data9/home/yzhao/GWAS_IRGSP1.0
stage=$1; pop=$2; folder=$3; filename=$4; runid=$5
[[ $pop =~ ^(GPall|GPallInd|GPallJap)$ && $folder =~ ^[A-Za-z0-9][A-Za-z0-9_.-]{0,79}$ ]] || exit 64
work="$BASE/$folder"
staging="$BASE/.gwas-web/$runid"
mkdir -p "$staging"
exec >> "$staging/run.log" 2>&1
PS4='+ $(date -Is) ${LINENO}: '
set -x
finish() {
  rc=$?
  trap - EXIT
  if ((rc != 0)) && [[ -d $work ]]; then
    find "$work" -maxdepth 1 -type f -name '*.sh.*' -exec tail -n 100 {} \; || true
  fi
  printf '%s\n' "$rc" > "$staging/$stage.exit.tmp"
  mv "$staging/$stage.exit.tmp" "$staging/$stage.exit"
}
trap finish EXIT
die() { echo "ERROR: $*"; exit 1; }
command -v perl >/dev/null
[[ $(id -un) == yzhao ]] || die 'Expected remote user yzhao'
stem=${filename%.*}
src="$BASE/000data_prepare/EMMAx.Data/${pop}_miss20"
if [[ $stage == prepare ]]; then
  [[ ! -e $work ]] || die "Directory already exists: $work (will not overwrite or resubmit)"
  [[ -s $src/$pop.sampleList ]] || die "Missing $src/$pop.sampleList"
  [[ -s $src/${pop}_allChr_snpNonHet.hIBS.kinf ]] || die 'Missing kinship file'
  [[ -d $BASE/000data_prepare/program ]] || die 'Missing program directory'
  shopt -s nullglob
  geno=("$src/$pop"*.tped)
  ((${#geno[@]})) || die 'No genotype tped files'
  mkdir "$work"
  cd "$work"
  mkdir chromosome EMMAx.Data EMMAx.Result phenotype plot QQplot sh.file sigSNP tped
  cp -r "$BASE/000data_prepare/program" .
  cp "$src/$pop.sampleList" phenotype/
  cp "${geno[@]}" tped/
  tfams=("$src/$pop"*.tfam)
  if ((${#tfams[@]})); then cp "${tfams[@]}" tped/; fi
  cp "$src/${pop}_allChr_snpNonHet.hIBS.kinf" tped/
  cp "$staging/phenotype.upload" "phenotype/$filename"
  for tool in select.pheno.only.pl EMMAx-Step2.ChangetoEMMAxPhenotype.tfam.pl qsub_nodeAdmin.pl EMMAx.run-get.pValue.hIBS.pl plot_Allpicture_ofOneTrait.pl stat.allSNPSignificantPos.inOneDir.list.pl; do
    [[ -s program/$tool ]] || die "Missing program/$tool"
  done
  perl program/select.pheno.only.pl "phenotype/$pop.sampleList" "phenotype/$filename" "phenotype/${pop}_${stem}.csv"
  [[ -s phenotype/${pop}_${stem}.csv ]] || die 'Empty extracted phenotype'
  perl program/EMMAx-Step2.ChangetoEMMAxPhenotype.tfam.pl "phenotype/${pop}_${stem}.csv" "$pop"
elif [[ $stage == admin ]]; then
  cd "$work"
  command -v qstat >/dev/null
  command -v qsub >/dev/null
  mkdir -p .gwas-status
  # Each queued command records its actual exit status on the shared filesystem.
  cat > .gwas-status/worker.sh <<'WORKER'
#!/usr/bin/env bash
set -Eeuo pipefail
cd "$1"
key=$2
shift 2
done_marker() {
  rc=$?
  trap - EXIT
  printf '%s\n' "$rc" > ".gwas-status/$key.exit.tmp"
  mv ".gwas-status/$key.exit.tmp" ".gwas-status/$key.exit"
}
trap done_marker EXIT
"$@"
WORKER
  submit() {
    local key=$1 name=$2 command output id
    shift 2
    printf -v command '%q ' bash "$work/.gwas-status/worker.sh" "$work" "$key" "$@"
    output=$(perl program/qsub_nodeAdmin.pl "$name" "$command") || die "Submission failed: $name"
    printf '%s\n' "$output"
    # PBS/Torque: 123.server ; SGE: Your job 123 (...) has been submitted.
    id=$(printf '%s\n' "$output" | sed -nE 's/^([0-9]+)(\.[A-Za-z0-9_.-]+)?[[:space:]]*$/\1/p; s/.*Your job(-array)? ([0-9]+).*/\2/p' | sort -u)
    [[ $id =~ ^[0-9]+$ ]] || die "Cannot identify exactly one queue job ID for $name. Job may have been submitted; inspect queue before retrying. Output: $output"
    printf '%s\t%s\n' "$key" "$id" >> ".gwas-status/$phase.jobs"
  }
  wait_jobs() {
    local start=$SECONDS snapshot key id code all active missing_since=0
    while true; do
      snapshot=$(qstat -u yzhao) || die 'qstat failed; refusing to infer completion'
      printf '%s\n' "$snapshot" > .gwas-status/queue.latest
      all=1; active=0
      while read -r key id; do
        if printf '%s\n' "$snapshot" | awk -v id="$id" '{split($1,a,"."); if(a[1]==id) found=1} END {exit !found}'; then active=1; fi
        if [[ -f .gwas-status/$key.exit ]]; then
          code=$(cat ".gwas-status/$key.exit")
          [[ $code == 0 ]] || die "Job $id ($key) failed with exit code $code; inspect *.sh.* logs"
        else all=0; fi
      done < ".gwas-status/$phase.jobs"
      if ((all == 1 && active == 0)); then break; fi
      if ((all == 0 && active == 0)); then
        if ((missing_since == 0)); then missing_since=$SECONDS; fi
        ((SECONDS-missing_since < 180)) || die 'Jobs disappeared without completion markers (killed, scheduler error, or unshared filesystem)'
      else missing_since=0; fi
      ((SECONDS-start < 604800)) || die 'Monitoring exceeded 7 days; jobs are not cancelled'
      echo "WAIT: $phase; completion_markers=$all; jobs_in_queue=$active; elapsed=$((SECONDS-start))s"
      sleep 20
    done
    # Scheduler output files may appear shortly after jobs leave the queue.
    sleep 10
    shopt -s nullglob
    logs=(*.sh.*)
    ((${#logs[@]})) || die 'No scheduler logs found; cannot verify errors'
    if grep -inE 'error|fatal|segmentation fault|traceback' "${logs[@]}"; then
      die "Error text detected after $phase; logs preserved in working directory"
    else
      code=$?
      [[ $code == 1 ]] || die 'Failed to read scheduler logs'
    fi
  }
  shopt -s nullglob
  phase=emmax
  files=(EMMAx.Data/"${pop}_"*.tped)
  ((${#files[@]})) || die 'No EMMAx.Data tped files generated by phenotype conversion'
  # Validate every input before the first submission.
  for f in "${files[@]}"; do
    name=${f##*/}; name=${name%.tped}
    [[ $name =~ ^(${pop}_[A-Za-z0-9_.-]+)_chr([0-9]+)$ ]] || die "Unexpected tped name: $name"
    trait=${BASH_REMATCH[1]}; chr=${BASH_REMATCH[2]}
    [[ -s EMMAx.Data/$trait.pheno || -s $trait.pheno ]] || die "Missing phenotype $trait.pheno"
    [[ -s /public/home/yzhao/IRGSP1.0_GPall/MAF/${pop}_chr${chr}_MAF.vcf ]] || die "Missing MAF for chromosome $chr"
  done
  : > .gwas-status/emmax.jobs
  for f in "${files[@]}"; do
    name=${f##*/}; name=${name%.tped}; trait=${name%_chr*}; chr=${name##*_chr}
    submit "${name}_emmax" "${name}_emmax_hibs.sh" perl program/EMMAx.run-get.pValue.hIBS.pl "$name" "$trait.pheno" "tped/${pop}_allChr_snpNonHet.hIBS.kinf" "/public/home/yzhao/IRGSP1.0_GPall/MAF/${pop}_chr${chr}_MAF.vcf"
  done
  wait_jobs
  for f in "${files[@]}"; do
    name=${f##*/}; name=${name%.tped}
    [[ -s EMMAx.Result/hIBS/$name.hIBS ]] || die "Missing GWAS result $name.hIBS"
  done
  mkdir -p sh.file/emmax
  logs=(*_emmax_hibs.sh.o* *_emmax_hibs.sh.e*)
  ((${#logs[@]})) && mv "${logs[@]}" sh.file/emmax/
  phase=plot
  files=(EMMAx.Result/hIBS/"${pop}_"*_chr01.hIBS)
  ((${#files[@]})) || die 'No chr01 results for plotting'
  : > .gwas-status/plot.jobs
  for f in "${files[@]}"; do
    name=${f##*/}; trait=${name%_chr01.hIBS}
    submit "${trait}_plot" "${trait}_plotGWAS.sh" perl program/plot_Allpicture_ofOneTrait.pl "$f" "$trait" "$trait"
  done
  wait_jobs
  [[ -n $(find plot QQplot -type f -size +0c -print -quit) ]] || die 'No nonempty plot outputs'
  mkdir -p sh.file/plot
  logs=(*_plotGWAS.sh.o* *_plotGWAS.sh.e*)
  ((${#logs[@]})) && mv "${logs[@]}" sh.file/plot/
elif [[ $stage == final ]]; then
  cd "$work"
  perl program/stat.allSNPSignificantPos.inOneDir.list.pl EMMAx.Result/hIBS hIBS "$pop" "$pop"
  [[ -f sigSNP/$pop.list ]] || die 'Missing peak SNP list'
  perl /data5/home/yzhao/program/compared_peakSNP_with_RiceNaviGene_1.0.pl "sigSNP/$pop.list" "sigSNP/${pop}_geneList.csv"
  [[ -f sigSNP/${pop}_geneList.csv ]] || die 'Missing gene annotation output'
  echo 'SUCCESS: GWAS analysis completed'
else die 'Unknown stage'; fi
