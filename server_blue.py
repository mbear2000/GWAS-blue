#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GWAS-submit BLUE extension v4

Usage:
    Keep the original server.py and workflow_v2.sh in the same directory.
    Start this file instead of server.py.

Single phenotype:
    Uses the original server.py validation and workflow unchanged.

>=2 phenotype files:
    1. Clean each environment's trait headers.
    2. Infer Location + Year.
    3. Automatically select exactly one BLUE model/program:
       - one location / multiple years -> server generate_blue_multi_env_allTrait.R
       - multiple locations / one year -> server generate_blue_multi_env_allTrait.R
       - multiple locations / multiple years -> server generate_blue_multiEnvYear_allTrait.R
    4. The user-entered "本次标签" is preserved and passed to the R script
       as its second argument. The program never overwrites that label.
    5. workflow_blue_v8.sh computes that BLUE result.
    5. Each plan retains its original .txt and sample-filtered .csv.
    6. A combined GWAS phenotype is created with columns:
           <plan-label>_<trait>
       so the downstream EMMAx names become:
           <population>_<plan-label>_<trait>
    7. The remote copy of workflow_v2.sh is patched only for BLUE mode so it
       does NOT prepend one extra global batch label.

GPall / GPallInd / GPallJap use exactly the same logic; only their selected
data-source paths differ.
"""

import base64
import csv
import json
import re
import secrets
import subprocess
import time
from pathlib import Path

import server as core

ROOT = Path(__file__).resolve().parent
BLUE_WORKFLOW = ROOT / "workflow_blue_v8.sh"
ORIGINAL_WORKFLOW = ROOT / "workflow_v2.sh"

ORIGINAL_VALIDATE = core.validate
ORIGINAL_BATCH_LABEL = core.batch_label
ORIGINAL_PHENOTYPE_TRAITS = core.phenotype_traits
ORIGINAL_BRIDGE_SCRIPT = core.bridge_script
ORIGINAL_SNAPSHOT = core.snapshot

BUNDLE_MAGIC = b"#GWAS_BLUE_BUNDLE_V4\n"
MAX_FILE = 20 * 1024 * 1024


def _b64s(value: str) -> str:
    return base64.b64encode(value.encode("utf-8")).decode("ascii")


def _ub64(value: str) -> str:
    return base64.b64decode(value).decode("utf-8")


def _decode_header(content: bytes):
    try:
        first = content.decode("utf-8-sig").splitlines()[0]
        sep = "\t" if "\t" in first else ","
        cols = next(csv.reader([first], delimiter=sep))
    except (UnicodeError, IndexError, csv.Error) as exc:
        raise ValueError("表型须为 UTF-8 的 CSV/TSV/TXT 文本") from exc
    cols = [x.strip() for x in cols]
    if len(cols) < 2 or cols[0] != "Accession":
        raise ValueError("第一列表头必须为 Accession，且至少包含一个性状列")
    return cols


def _apply_cleanup(name: str, prefix: str) -> str:
    if prefix and name.startswith(prefix):
        name = name[len(prefix):]
    return name


def _clean_traits(content: bytes, prefix: str):
    cols = _decode_header(content)
    cleaned = [_apply_cleanup(x, prefix) for x in cols[1:]]

    for trait in cleaned:
        if not trait:
            raise ValueError("清理后出现空性状名")
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,159}", trait):
            raise ValueError("清理后的性状名含不支持字符：" + trait)
        if "chr01" in trait:
            raise ValueError("性状名不能包含保留字 chr01：" + trait)

    if len(set(cleaned)) != len(cleaned):
        raise ValueError("清理后出现重复性状名，请检查前缀/额外删除文本")
    return cleaned


def _decode_file_item(item):
    if not isinstance(item, dict):
        raise ValueError("表型文件参数无效")

    filename = str(item.get("filename", "")).strip()
    prefix = str(item.get("prefix", "")).strip()

    if not re.fullmatch(
        r"[A-Za-z0-9][A-Za-z0-9_.-]{0,159}\.(txt|csv|tsv)", filename, re.I
    ):
        raise ValueError(
            "表型文件名请使用英文、数字、点、下划线或短横线，扩展名 txt/csv/tsv"
        )
    if prefix and not re.fullmatch(r"[A-Za-z0-9_.-]{1,140}", prefix):
        raise ValueError("环境前缀仅可含英文、数字、点、下划线或短横线")

    try:
        content = base64.b64decode(item.get("content", ""), validate=True)
    except Exception as exc:
        raise ValueError("表型文件编码无效：" + filename) from exc

    if not content or len(content) > MAX_FILE or b"\0" in content:
        raise ValueError("表型文件为空、过大或不是文本：" + filename)

    content = (
        content.removeprefix(b"\xef\xbb\xbf")
        .replace(b"\r\n", b"\n")
        .replace(b"\r", b"\n")
    )

    return {
        "filename": filename,
        "prefix": prefix,
        "content": content,
    }


def _infer_environment(item):
    """
    Infer location and year from:
      1) environment prefix
      2) first trait header
      3) filename

    Supports examples:
      LSh2022_549_
      SJ2023_
      GPallInd_LSh2023_Days-Height_...
    """
    cols = _decode_header(item["content"])
    candidates = [
        item["prefix"].rstrip("_"),
        cols[1],
        Path(item["filename"]).stem,
    ]

    patterns = [
        re.compile(r"^(?:GPallInd_|GPallJap_|GPall_)?([A-Za-z]+)(20\d{2})(?:_|$)"),
        re.compile(r"(?:^|_)([A-Za-z]+)(20\d{2})(?:_|$)"),
    ]

    for source in candidates:
        for pat in patterns:
            m = pat.search(source)
            if m:
                return m.group(1), m.group(2)

    raise ValueError(
        f"无法从 {item['filename']} 识别地点和年份。"
        "请让文件名/环境前缀/表头包含类似 LSh2022 或 SJ2024。"
    )


def _prepare_items(data):
    phenotypes = data.get("phenotypes")
    if not isinstance(phenotypes, list) or len(phenotypes) < 2:
        raise ValueError("BLUE 模式至少需要2个表型文件")

    items = [_decode_file_item(x) for x in phenotypes]
    names = [x["filename"] for x in items]
    if len(names) != len(set(names)):
        raise ValueError("不能重复添加同名表型文件")

    loc_case = {}
    seen_env = set()
    baseline_traits = None

    for item in items:
        loc, year = _infer_environment(item)
        loc = loc_case.setdefault(loc.lower(), loc)
        item["location"] = loc
        item["year"] = year
        item["traits"] = _clean_traits(
            item["content"], item["prefix"]
        )

        if baseline_traits is None:
            baseline_traits = item["traits"]
        elif item["traits"] != baseline_traits:
            raise ValueError(
                "多个环境清理后的性状名称或顺序不一致。\n"
                f"基准：{', '.join(baseline_traits[:8])}\n"
                f"{item['filename']}：{', '.join(item['traits'][:8])}"
            )

        env_key = (loc.lower(), year)
        if env_key in seen_env:
            raise ValueError(
                f"重复环境：{loc}{year}。每个地点×年份只能上传一个表型文件。"
            )
        seen_env.add(env_key)

    return items, baseline_traits


def _analysis_stem(traits):
    # Historical naming for a single corrected/special trait includes the
    # trait in the BLUE output filename. Multi-trait BLUE files omit it.
    return traits[0] if len(traits) == 1 else ""



def _blue_user_label(data):
    label = str(data.get("label", "")).strip()
    if not label:
        raise ValueError("请填写“本次标签”")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,79}", label):
        raise ValueError(
            "本次标签仅可包含英文、数字、点、下划线和短横线，最多80字符"
        )
    return label


def _filename_label(label):
    # Historical BLUE output filenames use underscores where analysis labels
    # commonly use hyphens, e.g. LocYear-blue -> LocYear_blue.
    return label.replace("-", "_")


def _automatic_plan(items, pop, traits, user_label):
    """
    Select exactly one BLUE model/program from environment structure.
    The analysis label is always supplied by the user.

    1 location + >=2 years:
        multi_env program

    >=2 locations + 1 year:
        multi_env program

    >=2 locations + >=2 years:
        locyear program

    user_label is passed unchanged to the R script. Only the output filename
    replaces "-" with "_" for compatibility with the user's historical naming.
    """
    locations = sorted({x["location"] for x in items}, key=str.lower)
    years = sorted({x["year"] for x in items})

    stem = _analysis_stem(traits)
    stem_part = (stem + "_") if stem else ""
    filename_label = _filename_label(user_label)

    if len(locations) == 1 and len(years) >= 2:
        return {
            "id": "multi_env_single_location",
            "kind": "multi_env",
            "label": user_label,
            "output": f"{pop}_{stem_part}{filename_label}",
            "title": f"{locations[0]} 多年单点 BLUE",
            "indexes": list(range(len(items))),
        }

    if len(locations) >= 2 and len(years) == 1:
        return {
            "id": "multi_env_single_year",
            "kind": "multi_env",
            "label": user_label,
            "output": f"{pop}_{stem_part}{filename_label}",
            "title": f"{years[0]} 单年多点 BLUE",
            "indexes": list(range(len(items))),
        }

    if len(locations) >= 2 and len(years) >= 2:
        if len(items) < 3:
            raise ValueError(
                "检测到多个地点和多个年份，但只有2个环境文件；"
                "Location 与 Year 会完全混杂，不能进行多年多点 BLUE。"
            )
        return {
            "id": "locyear",
            "kind": "locyear",
            "label": user_label,
            "output": f"{pop}_{stem_part}{filename_label}",
            "title": "多年多点 Loc × Year BLUE",
            "indexes": list(range(len(items))),
        }

    raise ValueError(
        "当前表型文件不能形成 BLUE。"
        "至少需要同地点多个年份、同年份多个地点，或多年多点数据。"
    )


def _folder_and_log(data, pop, filename_for_log):
    folder_input = str(data.get("directory", "")).strip()
    folder_suffix = re.sub(
        r"^(?:GPallInd|GPallJap|GPall)?miss20_", "", folder_input
    )
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,79}", folder_suffix):
        raise ValueError(
            "远程目录名称只填写 miss20_ 后面的目录主体；"
            "仅含字母、数字、下划线、点和短横线"
        )
    folder = f"{pop}miss20_{folder_suffix}"
    if len(folder) > 100:
        raise ValueError("远程目录名称过长")

    log_dir = Path(str(data.get("logDirectory", "")))
    if not log_dir.is_absolute() or not log_dir.is_dir():
        raise ValueError("本地日志目录必须是已经存在的绝对路径")

    log = log_dir / (Path(filename_for_log).stem + "_GWASrun_" + pop + ".log")
    return folder, log


def _encode_bundle(items, plans, pop, traits):
    rows = [BUNDLE_MAGIC.rstrip(b"\n")]
    rows.append(b"META\tpopulation\t" + _b64s(pop).encode("ascii"))

    for i, trait in enumerate(traits):
        rows.append(
            b"TRAIT\t%d\t%s" % (i, _b64s(trait).encode("ascii"))
        )

    for i, item in enumerate(items):
        rows.append(
            b"FILE\t%d\t%s\t%s\t%s\t%s"
            % (
                i,
                _b64s(item["filename"]).encode("ascii"),
                _b64s(item["prefix"]).encode("ascii"),
                _b64s(item["location"]).encode("ascii"),
                _b64s(item["year"]).encode("ascii"),
            )
        )
        encoded = base64.b64encode(item["content"]).decode("ascii")
        for start in range(0, len(encoded), 3000):
            rows.append(
                f"DATA\t{i}\t{encoded[start:start+3000]}".encode("ascii")
            )
        rows.append(f"END\t{i}".encode("ascii"))

    for i, plan in enumerate(plans):
        indexes = ",".join(str(x) for x in plan["indexes"])
        rows.append(
            b"PLAN\t%d\t%s\t%s\t%s\t%s\t%s\t%s"
            % (
                i,
                _b64s(plan["id"]).encode("ascii"),
                _b64s(plan["kind"]).encode("ascii"),
                _b64s(plan["label"]).encode("ascii"),
                _b64s(plan["output"]).encode("ascii"),
                _b64s(plan["title"]).encode("ascii"),
                _b64s(indexes).encode("ascii"),
            )
        )

    return b"\n".join(rows) + b"\n"


def _parse_bundle(content: bytes):
    if not content.startswith(BUNDLE_MAGIC):
        return None

    meta, files, chunks, traits, plans = {}, {}, {}, [], []

    for raw in content.decode("ascii").splitlines()[1:]:
        p = raw.split("\t")
        kind = p[0] if p else ""

        if kind == "META" and len(p) == 3:
            meta[p[1]] = _ub64(p[2])

        elif kind == "TRAIT" and len(p) == 3:
            traits.append(_ub64(p[2]))

        elif kind == "FILE" and len(p) == 6:
            idx = int(p[1])
            files[idx] = {
                "filename": _ub64(p[2]),
                "prefix": _ub64(p[3]),
                "location": _ub64(p[4]),
                "year": _ub64(p[5]),
            }
            chunks[idx] = []

        elif kind == "DATA" and len(p) == 3:
            chunks[int(p[1])].append(p[2])

        elif kind == "END" and len(p) == 2:
            pass

        elif kind == "PLAN" and len(p) == 8:
            plans.append({
                "id": _ub64(p[2]),
                "kind": _ub64(p[3]),
                "label": _ub64(p[4]),
                "output": _ub64(p[5]),
                "title": _ub64(p[6]),
                "indexes": [
                    int(x) for x in _ub64(p[7]).split(",") if x != ""
                ],
            })

        else:
            raise ValueError("BLUE bundle 格式无效")

    ordered = []
    for idx in sorted(files):
        item = files[idx]
        item["content"] = base64.b64decode("".join(chunks[idx]))
        ordered.append(item)

    if len(ordered) < 2 or not traits or not plans:
        raise ValueError("BLUE bundle 缺少必要内容")

    return {
        "meta": meta,
        "files": ordered,
        "traits": traits,
        "plans": plans,
    }


def phenotype_traits(content, pop, label):
    parsed = _parse_bundle(content)
    if parsed is None:
        return ORIGINAL_PHENOTYPE_TRAITS(content, pop, label)

    # Combined GWAS columns are <plan-label>_<clean-trait>.
    result = []
    for plan in parsed["plans"]:
        for trait in parsed["traits"]:
            result.append(f"{pop}_{plan['label']}_{trait}")

    if len(result) != len(set(result)):
        raise ValueError("BLUE 组合产生了重复 GWAS 性状名")
    return result


def batch_label(data, filename):
    phenotypes = data.get("phenotypes")
    if isinstance(phenotypes, list) and len(phenotypes) >= 2:
        return _blue_user_label(data)
    return ORIGINAL_BATCH_LABEL(data, filename)


def validate(data):
    phenotypes = data.get("phenotypes")

    # Single phenotype = unchanged upstream behavior.
    if not isinstance(phenotypes, list) or len(phenotypes) <= 1:
        if isinstance(phenotypes, list) and len(phenotypes) == 1:
            one = phenotypes[0]
            translated = dict(data)
            translated["filename"] = one.get("filename", "")
            translated["content"] = one.get("content", "")
            translated.pop("phenotypes", None)
            return ORIGINAL_VALIDATE(translated)
        return ORIGINAL_VALIDATE(data)

    pop = str(data.get("population", ""))
    if pop not in ("GPall", "GPallInd", "GPallJap"):
        raise ValueError("请选择有效群体")

    user_label = _blue_user_label(data)
    items, traits = _prepare_items(data)
    plans = [_automatic_plan(items, pop, traits, user_label)]

    # A synthetic filename is used only as the one combined phenotype passed
    # into the existing workflow. Individual BLUE files keep historical names.
    filename = f"{pop}_{_filename_label(user_label)}.txt"
    folder, log = _folder_and_log(data, pop, filename)
    bundle = _encode_bundle(items, plans, pop, traits)

    # Run the same trait validation that preview/run creation relies on.
    phenotype_traits(bundle, pop, "BLUE")
    return pop, folder, filename, bundle, log


def _upload_bytes_commands(remote_path, payload):
    q = core.shlex.quote
    encoded = base64.b64encode(payload.replace(b"\r\n", b"\n")).decode("ascii")
    cmds = [": > " + q(remote_path + ".b64")]
    for i in range(0, len(encoded), 3000):
        cmds.append(
            "printf '%s' '"
            + encoded[i:i+3000]
            + "' >> "
            + q(remote_path + ".b64")
        )
    cmds += [
        "base64 -d " + q(remote_path + ".b64") + " > " + q(remote_path),
        "test -s " + q(remote_path),
    ]
    return cmds


def _patched_original_workflow():
    """
    Patch the real upstream workflow_v2.sh only for BLUE runs.

    BLUE naming:
      filename/archive stem:
          <population>_<user-label-with-hyphen-as-underscore>

      final GWAS trait:
          <population>_<original-user-label>_<trait>

    Non-BLUE behavior is unchanged.
    """
    text = ORIGINAL_WORKFLOW.read_text(encoding="utf-8")

    trait_pattern = re.compile(
        r'(?P<indent>[ \t]*)if \[\[ \$old == "\$\{pop\}_\$\{tag\}" \|\| '
        r'\$old == "\$\{pop\}_\$\{tag\}_"\* \]\]; then trait=\$old\s*\n'
        r'(?P=indent)else trait="\$\{pop\}_\$\{tag\}_\$\{old#\$\{pop\}_\}"; fi'
    )
    trait_replacement = (
        r'\g<indent>if [[ -s "$staging/blue-mode.flag" ]]; then\n'
        r'\g<indent>  trait=$old\n'
        r'\g<indent>elif [[ $old == "${pop}_${tag}" || $old == "${pop}_${tag}_"* ]]; then\n'
        r'\g<indent>  trait=$old\n'
        r'\g<indent>else\n'
        r'\g<indent>  trait="${pop}_${tag}_${old#${pop}_}"\n'
        r'\g<indent>fi'
    )
    text, n = trait_pattern.subn(trait_replacement, text, count=1)
    if n != 1:
        raise RuntimeError("无法定位 workflow_v2.sh 的 trait 重命名代码")

    archive_old = (
        'archive="${stem}__${runid}"\n'
        'if [[ -s $staging/archive-name.txt ]]; then\n'
        '  archive=$(cat "$staging/archive-name.txt")\n'
        '  [[ $archive =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*_[0-9]{8}-[0-9]{4}$ ]] || exit 64\n'
        'fi'
    )
    archive_new = (
        'if [[ -s "$staging/blue-mode.flag" ]]; then\n'
        '  archive="$stem"\n'
        'else\n'
        '  archive="${stem}__${runid}"\n'
        '  if [[ -s $staging/archive-name.txt ]]; then\n'
        '    archive=$(cat "$staging/archive-name.txt")\n'
        '    [[ $archive =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*_[0-9]{8}-[0-9]{4}$ ]] || exit 64\n'
        '  fi\n'
        'fi'
    )
    if archive_old not in text:
        raise RuntimeError("无法定位 workflow_v2.sh 的 archive 初始化代码")
    text = text.replace(archive_old, archive_new, 1)

    collision_pattern = re.compile(
        r"(?m)^  \[\[ ! -d sh\.file/emmax/\$archive && ! -d sh\.file/plotGWAS/\$archive \]\] "
        r"\|\| die 'Queue-log archive name already exists; choose a new run time'\s*$"
    )
    collision_new = (
        '  if [[ -s "$staging/blue-mode.flag" ]]; then\n'
        '    rm -rf "sh.file/emmax/$archive" "sh.file/plotGWAS/$archive"\n'
        '  else\n'
        "    [[ ! -d sh.file/emmax/$archive && ! -d sh.file/plotGWAS/$archive ]] || "
        "die 'Queue-log archive name already exists; choose a new run time'\n"
        '  fi'
    )
    text, n = collision_pattern.subn(collision_new, text, count=1)
    if n != 1:
        raise RuntimeError("无法定位 workflow_v2.sh 的 archive 冲突检查代码")

    peak_old = '  peak_prefix="${pop}_${stem}"'
    peak_new = (
        '  if [[ -s "$staging/blue-mode.flag" ]]; then\n'
        '    peak_prefix="$stem"\n'
        '  else\n'
        '    peak_prefix="${pop}_${stem}"\n'
        '  fi'
    )
    if peak_old not in text:
        raise RuntimeError("无法定位 workflow_v2.sh 的 peak_prefix 代码")
    text = text.replace(peak_old, peak_new, 1)

    dest_old = (
        '  dest="$work/EMMAx.Result/hIBS/$stem"\n'
        '  [[ ! -e $dest || -d $dest ]] || die "Refusing to replace non-directory result archive: $dest"\n'
        '  mkdir -p "$dest/output"'
    )
    dest_new = (
        '  dest="$work/EMMAx.Result/hIBS/$stem"\n'
        '  if [[ -s "$staging/blue-mode.flag" && -d "$dest" ]]; then\n'
        '    rm -rf "$dest"\n'
        '  fi\n'
        '  [[ ! -e $dest || -d $dest ]] || die "Refusing to replace non-directory result archive: $dest"\n'
        '  mkdir -p "$dest/output"'
    )
    if dest_old not in text:
        raise RuntimeError("无法定位 workflow_v2.sh 的 result archive 代码")
    text = text.replace(dest_old, dest_new, 1)

    return text.encode("utf-8")

def bridge_script(run, pop, folder, filename, label=None, source=None):
    script = ORIGINAL_BRIDGE_SCRIPT(
        run, pop, folder, filename, label, source
    )
    uploaded = (run / "phenotype.upload").read_bytes()
    parsed = _parse_bundle(uploaded)

    q = core.shlex.quote
    work = core.BASE + "/" + folder
    remote = core.BASE + "/.gwas-web/" + run.name

    # Single-phenotype mode still runs through workflow_blue_v8.sh.
    # That wrapper delegates one-file analysis to a sibling workflow_v2.sh,
    # so the original workflow_v2.sh must also be uploaded for single-file
    # submissions. Previous versions only uploaded it for BLUE bundles.
    if parsed is None:
        if not ORIGINAL_WORKFLOW.is_file():
            raise RuntimeError("本地缺少 workflow_v2.sh")

        blocks = [
            "  Call Status(" + core.vb("Uploading original workflow_v2.sh") + ")"
        ]
        for command in _upload_bytes_commands(
            remote + "/workflow_v2.sh",
            ORIGINAL_WORKFLOW.read_bytes(),
        ):
            blocks.append(
                "  Call ExecChecked(" + core.vb(command) + ")"
            )

        marker = '  Call Status("Starting prepare stage")'
        if marker not in script:
            raise RuntimeError("无法定位 prepare 启动位置")

        script = script.replace(
            marker,
            "\n".join(blocks) + "\n" + marker,
            1,
        )
        return script

    # One overwrite confirmation for all uploaded inputs and all selected
    # BLUE outputs.
    conflicts = [
        work + "/phenotype/" + item["filename"]
        for item in parsed["files"]
    ]
    for plan in parsed["plans"]:
        conflicts += [
            work + "/lme4/blue/" + plan["output"] + ".txt",
            work + "/lme4/blue/" + plan["output"] + ".csv",
        ]
    conflicts.append(work + "/phenotype/" + filename)

    archive_stem = Path(filename).stem
    conflicts += [
        work + "/sh.file/emmax/" + archive_stem,
        work + "/sh.file/plotGWAS/" + archive_stem,
        work + "/EMMAx.Result/hIBS/" + archive_stem,
        work + "/sigSNP/" + archive_stem + ".list",
        work + "/sigSNP/" + archive_stem + "_geneList.csv",
    ]

    # Do not use "exit" in this scalar command. ExecScalar still needs the
    # current shell to print its scalar/end markers after the check returns.
    tests = " || ".join(f"test -e {q(path)}" for path in conflicts)
    checks = f"if {tests}; then printf yes; else printf no; fi"

    replacement = (
        "  duplicatePhenotype = ExecScalar(" + core.vb(checks) + ")"
    )
    script, n = re.subn(
        r'(?m)^  duplicatePhenotype = ExecScalar\(.+\)$',
        replacement,
        script,
        count=1,
    )
    if n != 1:
        raise RuntimeError("无法定位 SecureCRT 同名文件检查代码")

    names = "、".join(x["filename"] for x in parsed["files"])
    outputs = "、".join(x["output"] for x in parsed["plans"])
    msg = (
        f"BLUE 输入或输出已经存在。输入：{names}；输出：{outputs}。"
        "请在网页确认是否覆盖本次同名目标。"
    )
    script, n = re.subn(
        r'(?m)^    Call Save\("state\.txt", "confirm\|.*$',
        '    Call Save("state.txt", ' + core.vb("confirm|" + msg) + ")",
        script,
        count=1,
    )
    if n != 1:
        raise RuntimeError("无法定位 SecureCRT 覆盖提示代码")

    # Upload only the patched workflow_v2.sh.
    # BLUE R analysis scripts are no longer uploaded from Windows.
    # workflow_blue_v8.sh copies the fixed R scripts directly from:
    # /data9/home/yzhao/program/EMMAx/
    if not ORIGINAL_WORKFLOW.is_file():
        raise RuntimeError("本地缺少 workflow_v2.sh")

    payloads = [
        (
            remote + "/workflow_v2.sh",
            _patched_original_workflow(),
            "Uploading BLUE-aware workflow_v2.sh",
        ),
    ]

    blocks = []
    for remote_path, payload, status in payloads:
        blocks.append("  Call Status(" + core.vb(status) + ")")
        for command in _upload_bytes_commands(remote_path, payload):
            blocks.append(
                "  Call ExecChecked(" + core.vb(command) + ")"
            )

    marker = '  Call Status("Starting prepare stage")'
    if marker not in script:
        raise RuntimeError("无法定位 prepare 启动位置")

    script = script.replace(
        marker,
        "\n".join(blocks) + "\n" + marker,
        1,
    )
    return script



def _local_state_text(run):
    path = run / "state.txt"
    if not path.exists():
        return ""
    try:
        raw = path.read_bytes()
        enc = "utf-16" if raw.startswith((b"\\xff\\xfe", b"\\xfe\\xff")) else "utf-8-sig"
        return raw.decode(enc, errors="replace").strip()
    except OSError:
        return ""


def blue_snapshot(run):
    snap = ORIGINAL_SNAPSHOT(run)
    flag = run / "cancelled.flag"
    if flag.exists():
        try:
            reason = flag.read_text("utf-8").strip()
        except OSError:
            reason = ""
        snap["status"] = "cancelled"
        snap["detail"] = reason or "用户已停止本地网页监控；远程队列/进程未自动 qdel。"
    return snap


def _inflation_bridge_script(run, meta):
    pop = meta["population"]
    folder = meta["directory"]
    filename = meta["filename"]
    label = meta.get("label") or Path(filename).stem
    source = meta.get("dataSource")

    s = ORIGINAL_BRIDGE_SCRIPT(run, pop, folder, filename, label, source)

    q = core.shlex.quote
    remote = core.BASE + "/.gwas-web/" + run.name
    workflow_remote = remote + "/workflow_inflation.sh"

    upload_lines = []
    payload = ORIGINAL_WORKFLOW.read_bytes()
    for command in _upload_bytes_commands(workflow_remote, payload):
        upload_lines.append("  Call ExecChecked(" + core.vb(command) + ")")
    uploads = "\n".join(upload_lines)

    cmd = " ".join(
        q(x) for x in [
            "bash", workflow_remote, "inflation", pop, folder,
            filename, run.name, label
        ]
    )
    launch = "nohup " + cmd + " </dev/null >/dev/null 2>&1 &"

    main = (
        'Sub Main\n'
        '  Call Save("state.txt", "connecting|admin2 -> fat2；准备仅计算膨胀系数")\n'
        '  Set tab = crt.Session.ConnectInTab("/S admin2")\n'
        '  tab.Screen.Synchronous = True\n'
        '  tab.Screen.IgnoreEscape = True\n'
        '  Call PromptReady()\n'
        '  Call QuietPrompt("Connected to admin2")\n\n'
        '  tab.Screen.Send "ssh fat2" & vbCr\n'
        '  Call PromptReady()\n'
        '  Call QuietPrompt("Connected to fat2")\n'
        + '  Call ExecChecked("test $(id -un) = yzhao && test $(hostname -s) = fat2 && test -d ' + core.BASE + '")\n\n'
        '  Call Progress("Inflation-only recovery: existing archived .ps files will be reused")\n'
        + '  Call ExecChecked("mkdir -p ' + remote + '")\n'
        + uploads + '\n'
        + '  Call ExecChecked("rm -f ' + remote + '/inflation.exit")\n'
        + '  Call Stage("inflation", ' + core.vb(launch) + ')\n'
        '  Call Save("state.txt", "complete|膨胀系数计算完成；未重新运行 prepare/qsub/plot/GWAS")\n'
        '  Call Progress("Inflation factor calculation completed successfully")\n'
        '  tab.Screen.Send "stty echo; exit" & vbCr\n'
        'End Sub'
    )

    begin = s.index("Sub Main")
    end = s.index("\nEnd Sub", begin) + len("\nEnd Sub")
    return s[:begin] + main + s[end:]


class BlueHandler(core.Handler):
    def do_POST(self):
        if self.path not in ("/api/inflation-only", "/api/cancel-local"):
            return super().do_POST()

        if not self.valid_host() or self.headers.get("X-GWAS-Token") != core.TOKEN:
            return self.send({"error": "请求校验失败，请刷新页面"}, 403)

        try:
            length = int(self.headers.get("Content-Length", "0"))
            if not 0 < length <= 1024 * 1024:
                raise ValueError("请求过大或为空")
            data = json.loads(self.rfile.read(length))
            rid = str(data.get("id", ""))
            if not re.fullmatch(r"[A-Za-z0-9_-]+", rid):
                raise ValueError("无效任务编号")

            with core.LOCK:
                old_run = core.RUNS / rid
                if not (old_run / "meta.json").exists():
                    raise ValueError("本地任务不存在")
                old_meta = json.loads((old_run / "meta.json").read_text("utf-8"))

                if self.path == "/api/cancel-local":
                    reason = (
                        "用户已中断/清除本地网页任务状态。"
                        "此操作不会自动 qdel 服务器队列，也不会删除服务器结果。"
                    )
                    (old_run / "cancelled.flag").write_text(reason, "utf-8")
                    (old_run / "state.txt").write_text("cancelled|" + reason, "utf-8")
                    return self.send(blue_snapshot(old_run))

                if old_meta.get("remoteDirectoryDeleted"):
                    raise ValueError("服务器工作目录已删除，不能从已有 .ps 恢复膨胀系数")

                for prior in core.RUNS.iterdir():
                    if prior == old_run or not (prior / "meta.json").exists():
                        continue
                    ps = blue_snapshot(prior)
                    if ps["status"] in ("running", "connecting"):
                        raise ValueError("另一个网页任务正在运行。请先停止/核实该任务，再启动“只计算膨胀系数”。")

                if not core.CRT.exists():
                    raise ValueError("未找到 SecureCRT: " + str(core.CRT))

                old_reason = "已由“只计算膨胀系数”恢复任务接管；旧网页监控状态已停止。"
                (old_run / "cancelled.flag").write_text(old_reason, "utf-8")
                (old_run / "state.txt").write_text("cancelled|" + old_reason, "utf-8")

                new_rid = time.strftime("%Y%m%d-%H%M%S-") + secrets.token_hex(4)
                run = core.RUNS / new_rid
                run.mkdir()

                src_upload = old_run / "phenotype.upload"
                if src_upload.exists():
                    (run / "phenotype.upload").write_bytes(src_upload.read_bytes())
                else:
                    (run / "phenotype.upload").write_bytes(b"placeholder\n")

                display = Path(old_meta["filename"]).stem + "_inflationFactor"
                old_log = Path(old_meta["logPath"])
                log = core.available_log_path(
                    old_log.with_name(display + "_" + time.strftime("%Y%m%d-%H%M%S") + ".log")
                )
                header = "\n".join([
                    "GWAS inflation-only run " + new_rid,
                    "开始时间: " + time.strftime("%Y-%m-%d %H:%M:%S"),
                    "来源任务: " + rid,
                    "群体: " + old_meta["population"],
                    "远程目录: " + core.BASE + "/" + old_meta["directory"],
                    "归档结果: EMMAx.Result/hIBS/" + Path(old_meta["filename"]).stem,
                    "模式: 仅使用已有 *.ps 计算 genomic inflation factor；不重新运行 GWAS",
                ])
                with log.open("x", encoding="utf-8") as f:
                    f.write(header)

                meta = dict(old_meta)
                meta.update(
                    id=new_rid,
                    logPath=str(log),
                    header=header,
                    displayName=display,
                    inflationOnlyOf=rid,
                )
                (run / "meta.json").write_text(json.dumps(meta, ensure_ascii=False), "utf-8")
                (run / "bridge.vbs").write_text(_inflation_bridge_script(run, meta), "utf-16")
                (run / "state.txt").write_text("connecting|启动 SecureCRT；只计算膨胀系数", "utf-8")

                startup = subprocess.STARTUPINFO()
                startup.dwFlags |= subprocess.STARTF_USESHOWWINDOW
                startup.wShowWindow = 0
                try:
                    subprocess.Popen([str(core.CRT), "/SCRIPT", str(run / "bridge.vbs")], startupinfo=startup)
                except OSError as exc:
                    (run / "state.txt").write_text("failed|SecureCRT inflation-only launch failed", "utf-8")
                    raise ValueError("SecureCRT 启动失败: " + str(exc)) from exc

                return self.send(blue_snapshot(run), 201)

        except (ValueError, OSError, TypeError, json.JSONDecodeError) as exc:
            return self.send({"error": str(exc)}, 400)



# Patch imported original server runtime.
core.validate = validate
core.batch_label = batch_label
core.phenotype_traits = phenotype_traits
core.bridge_script = bridge_script
core.snapshot = blue_snapshot
core.WORKFLOW = BLUE_WORKFLOW


if __name__ == "__main__":
    core.RUNS.mkdir(exist_ok=True)
    core.threading.Thread(
        target=core.mirror_logs, daemon=True
    ).start()
    print(
        f"GWAS BLUE v9.9 console: http://127.0.0.1:{core.PORT}",
        flush=True,
    )
    core.ThreadingHTTPServer(
        ("127.0.0.1", core.PORT), BlueHandler
    ).serve_forever()
