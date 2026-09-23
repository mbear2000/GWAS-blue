"""Local-only GWAS console. Python standard library; no credentials stored."""
import base64
import csv
import json
import os
from pathlib import Path
import re
import secrets
import shlex
import subprocess
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse, parse_qs

ROOT = Path(__file__).resolve().parent
RUNS = ROOT / 'runs'
CRT = Path(r'D:\software\SecureCRT\App\VanDyke Clients\SecureCRT.exe')
PORT = int(os.environ.get('GWAS_WEB_PORT', '8765'))
TOKEN = secrets.token_urlsafe(32)
LOCK = threading.Lock()
NOTIFY_LOCK = threading.Lock()
BASE = '/data9/home/yzhao/GWAS_IRGSP1.0'
WORKFLOW = ROOT / 'workflow_v2.sh'

def data_source(data, pop):
    src = BASE + '/000data_prepare/EMMAx.Data/' + pop + '_miss20/'
    defaults = dict(version=pop + 'miss20-864lines', sample=src + pop + '.sampleList',
                    tped=src + pop + '_chr{chr}.tped',
                    kinship=src + pop + '_allChr_snpNonHet.hIBS.kinf',
                    maf='/public/home/yzhao/IRGSP1.0_GPall/MAF/' + pop + '_chr{chr}_MAF.vcf')
    supplied = data.get('dataSource', {})
    if not isinstance(supplied, dict):
        raise ValueError('数据源配置无效')

    # A previously saved blank version must not suppress the current default.
    # Other user-supplied source fields remain fully editable.
    result = {}
    for k, default in defaults.items():
        value = supplied.get(k, default)
        if k == 'version' and isinstance(value, str) and not value.strip():
            value = default
        result[k] = value
    for k, v in result.items():
        if not isinstance(v, str):
            raise ValueError('数据源字段须为文本')
        if k == 'version':
            valid = re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,63}', v)
        else:
            clean = v.replace('{chr}', '01')
            valid = re.fullmatch(r'/[A-Za-z0-9_./-]+', clean) and '..' not in clean.split('/')
            valid = valid and (v.count('{chr}') == (1 if k in ('tped', 'maf') else 0))
        if not valid:
            raise ValueError('数据源字段无效：' + k + '；路径使用服务器绝对路径，TPED/MAF 各包含一次 {chr}')
    return result


def batch_label(data, filename):
    default = re.sub(r'_Re-fitting$', '', Path(filename).stem, flags=re.I)
    label = data.get('label', '').strip() or default
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,63}', label):
        raise ValueError('本次标签限64字符，使用英文、数字、点、下划线或短横线')
    return label


def phenotype_traits(content, pop, label):
    try:
        header = content.decode('utf-8-sig').splitlines()[0]
        columns = next(csv.reader([header], delimiter='\t' if '\t' in header else ','))
    except (UnicodeError, IndexError, csv.Error) as exc:
        raise ValueError('表型须为 UTF-8 的制表符或逗号分隔文本') from exc
    if len(columns) < 2 or columns[0] != 'Accession':
        raise ValueError('第一列表头必须为 Accession，后面至少一个性状列')
    names = columns[1:]
    if len(set(names)) != len(names):
        raise ValueError('表型文件内有重复的性状列名，请先区分列名')
    for name in names:
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,139}', name) or 'chr01' in name:
            raise ValueError('性状列名仅可含英文、数字、点、下划线、短横线，且不能含保留字 chr01：' + name)
    traits = [pop + '_' + (n if n == label or n.startswith(label + '_') else label + '_' + n) for n in names]
    if len(set(traits)) != len(traits) or any(len(n) > 170 or 'chr01' in n for n in traits):
        raise ValueError('加标签后的性状名称重复、过长或含保留字 chr01，请修改标签或列名')
    return traits


def validate(data):
    pop = data.get('population', '')
    folder_input = str(data.get('directory', '')).strip()
    filename = data.get('filename', '')
    if pop not in ('GPall', 'GPallInd', 'GPallJap'):
        raise ValueError('请选择有效群体')

    # The web field stores only the reusable suffix AFTER "<population>miss20_".
    # Examples accepted:
    #   UAV_549lines_5years_Re-fitting
    #   miss20_UAV_549lines_5years_Re-fitting
    #   GPallJapmiss20_UAV_549lines_5years_Re-fitting
    # The real directory is always rebuilt from the CURRENT population:
    #   <population>miss20_<suffix>
    folder_suffix = re.sub(r'^(?:GPallInd|GPallJap|GPall)?miss20_', '', folder_input)
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,79}', folder_suffix):
        raise ValueError('远程目录名称只填写 miss20_ 后面的目录主体；须以字母或数字开头，仅含字母、数字、下划线、点和短横线')
    folder = f'{pop}miss20_{folder_suffix}'
    if len(folder) > 80:
        raise ValueError('群体名称、miss20_ 与远程目录名称拼接后最多80字符，请缩短目录名称')
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,119}\.(txt|csv|tsv)', filename, re.I):
        raise ValueError('表型文件名请使用英文、数字、点、下划线或短横线，扩展名为 txt/csv/tsv')
    try:
        content = base64.b64decode(data.get('content', ''), validate=True)
    except Exception as exc:
        raise ValueError('文件编码无效') from exc
    if not content or len(content) > 20 * 1024 * 1024 or b'\0' in content:
        raise ValueError('请上传非空文本表型文件，最大20MB')
    content = content.removeprefix(b'\xef\xbb\xbf').replace(b'\r\n', b'\n')
    phenotype_traits(content, pop, batch_label(data, filename))
    directory = Path(data.get('logDirectory', ''))
    if not directory.is_absolute() or not directory.is_dir():
        raise ValueError('本地日志目录必须是已存在的绝对路径')
    log = directory / (Path(filename).stem + '_GWASrun_' + pop + '.log')
    return pop, folder, filename, content, log


def vb(value):
    return '"' + str(value).replace('"', '""') + '"'


def available_log_path(path):
    if not path.exists():
        return path
    stamp = time.strftime('%Y%m%d-%H%M%S')
    candidate = path.with_name(path.stem + '_' + stamp + path.suffix)
    counter = 2
    while candidate.exists():
        candidate = path.with_name(path.stem + '_' + stamp + '_' + str(counter) + path.suffix)
        counter += 1
    return candidate


def bridge_script(run, pop, folder, filename, label=None, source=None):
    rid = run.name
    remote = f'{BASE}/.gwas-web/{rid}'
    q = shlex.quote
    workflow = WORKFLOW.read_bytes().replace(b'\r\n', b'\n')
    label = label or batch_label({}, filename)
    uploaded = (run / 'phenotype.upload').read_bytes()
    expected = ('\n'.join(phenotype_traits(uploaded, pop, label)) + '\n').encode('utf-8')
    upload_commands = []
    payloads = [('workflow.sh', workflow), ('phenotype.upload', uploaded), ('expected_traits.txt', expected)]
    source = source or data_source({}, pop)
    payloads.append(('data-source.txt', ('\n'.join(source[k] for k in ('version', 'sample', 'tped', 'kinship', 'maf')) + '\n').encode()))
    for name, payload in payloads:
        encoded = base64.b64encode(payload).decode()
        upload_commands.append(f': > {q(remote + "/" + name + ".b64")}')
        for i in range(0, len(encoded), 3000):
            upload_commands.append(f"printf '%s' '{encoded[i:i+3000]}' >> {q(remote + '/' + name + '.b64')}")
        upload_commands.append(f'base64 -d {q(remote + "/" + name + ".b64")} > {q(remote + "/" + name)}')
    uploads = '\n'.join('  Call ExecChecked(' + vb(c) + ')' for c in upload_commands)
    launch = {}
    for stage in ('prepare', 'admin', 'final'):
        cmd = ' '.join(q(x) for x in ['bash', remote + '/workflow.sh', stage, pop, folder, filename, rid, label or batch_label({}, filename)])
        launch[stage] = f'nohup {cmd} </dev/null >/dev/null 2>&1 &'
    return f'''# $language = "VBScript"
# $interface = "1.0"
Option Explicit
Dim tab, fso, localDir, remoteDir, beginTag, endTag, scalarBegin, scalarEnd, logOffset, logText
logOffset = 0
logText = ""
localDir = {vb(run)}
remoteDir = {vb(remote)}
beginTag = "@B{rid[-8:]}@"
endTag = "@E{rid[-8:]}@"
scalarBegin = "@SB{rid[-8:]}@"
scalarEnd = "@SE{rid[-8:]}@"
Set fso = CreateObject("Scripting.FileSystemObject")

Sub Save(name, value)
  Dim f, path, attempt, errNum, errDesc
  path = localDir & "\\" & name
  errNum = 0
  errDesc = ""

  ' Do not use DeleteFile + MoveFile here.  server.py reads these small
  ' status files every few seconds, and Windows can briefly deny deletion
  ' while another process has the file open.  Overwrite in place and retry
  ' short transient sharing violations instead.
  For attempt = 1 To 30
    On Error Resume Next
    Err.Clear

    Set f = fso.CreateTextFile(path, True, True)
    If Err.Number = 0 Then
      f.Write value
      If Err.Number = 0 Then
        f.Close
        If Err.Number = 0 Then
          On Error GoTo 0
          Exit Sub
        End If
      Else
        f.Close
      End If
    End If

    errNum = Err.Number
    errDesc = Err.Description
    Err.Clear
    On Error GoTo 0

    ' SecureCRT provides crt.Sleep; 30 x 100 ms gives Windows enough time
    ' to release a transient read handle without blocking the workflow long.
    crt.Sleep 100
  Next

  Err.Raise vbObjectError + 101, "GWAS local status write", _
    "Cannot save " & path & " after 30 retries. Last error " & CStr(errNum) & ": " & errDesc
End Sub

Function Enc(s)
  Dim i, result
  result = ""
  For i = 1 To Len(s)
    result = result & "\\" & Right("000" & Oct(Asc(Mid(s,i,1))),3)
  Next
  Enc = result
End Function

Sub Fail(message)
  Call Save("state.txt", "attention|" & message)
  Err.Raise vbObjectError + 100, "GWAS", message
End Sub

Sub PromptReady()
  Dim result
  result = tab.Screen.WaitForStrings(Array("$ ", "# ", "> "), 300)
  If result = 0 Then Call Fail("Login not confirmed. Open SecureCRT and inspect authentication. No further commands sent.")
End Sub

Sub QuietTerminal()
  tab.Screen.Send "stty -echo" & vbCr
  Call PromptReady()
End Sub

Sub Status(message)
  tab.Screen.Send "printf '\\n[GWAS] %s\\n' " & ShellQuote(message) & vbCr
  Call PromptReady()
End Sub

Function ExecRemote(command)
  Dim output, result
  Call Save("last-command.txt", Left(command, 600))
  tab.Screen.Send "printf '" & Enc(beginTag) & "'; " & command & "; printf '" & Enc(endTag) & "'; printf '\\r\\033[2K'" & vbCr
  result = tab.Screen.WaitForString(beginTag, 120)
  If Not result Then Call Fail("Remote command start timeout; inspect SecureCRT. Jobs may still be running.")
  output = tab.Screen.ReadString(endTag, 120)
  If tab.Screen.MatchIndex = 0 Then Call Fail("Remote command timeout; inspect SecureCRT. Jobs may still be running.")
  Call PromptReady()
  ExecRemote = output
End Function

Sub ExecChecked(command)
  Dim output
  output = ExecRemote("bash -c " & ShellQuote(command) & " && printf '" & Enc("CHECK_OK") & "'")
  If InStr(output, "CHECK_OK") = 0 Then Call Fail("Remote setup failed: " & output)
  Call Save("heartbeat.txt", CStr(Now))
End Sub

Function ShellQuote(s)
  ShellQuote = "'" & Replace(s, "'", "'" & Chr(34) & "'" & Chr(34) & "'") & "'"
End Function

Function CleanValue(s)
  s = Replace(s, vbCr, "")
  s = Replace(s, vbLf, "")
  CleanValue = Trim(s)
End Function

Function ExecScalar(command)
  Dim raw, p1, p2, value
  raw = ExecRemote("printf '" & Enc(scalarBegin) & "'; " & command & "; printf '" & Enc(scalarEnd) & "'")
  p1 = InStr(1, raw, scalarBegin, vbBinaryCompare)
  If p1 = 0 Then Call Fail("Cannot find scalar begin marker in remote response.")
  p1 = p1 + Len(scalarBegin)
  p2 = InStr(p1, raw, scalarEnd, vbBinaryCompare)
  If p2 = 0 Then Call Fail("Cannot find scalar end marker in remote response.")
  value = Mid(raw, p1, p2 - p1)
  ExecScalar = CleanValue(value)
End Function

Sub SyncLog(drain)
  Dim size, count, chunks, piece
  size = ExecScalar("if test -f " & remoteDir & "/run.log; then wc -c < " & remoteDir & "/run.log; else printf 0; fi")
  If Not IsNumeric(size) Then Call Fail("Cannot read remote log size after scalar parsing: " & Left(size, 200))
  size = CDbl(size)
  If size < logOffset Then Call Fail("Remote log was truncated; inspect before continuing.")
  chunks = 0
  Do While logOffset < size
    count = size - logOffset
    If count > 4096 Then count = 4096
    piece = ExecRemote("dd if=" & remoteDir & "/run.log bs=1 skip=" & CStr(logOffset) & " count=" & CStr(count) & " 2>/dev/null")
    logText = logText & piece
    logOffset = logOffset + count
    Call Save("remote.log", logText)
    Call Save("heartbeat.txt", CStr(Now))
    chunks = chunks + 1
    If Not drain And chunks >= 4 Then Exit Do
  Loop
End Sub

Sub Stage(name, command)
  Call Save("state.txt", "running|" & name)
  Call ExecChecked(command & " wait_pid=$!; kill -0 $wait_pid")
  Call WaitStage(name)
End Sub

Sub WaitStage(name)
  Dim output, code, polls
  Call Save("state.txt", "running|" & name)
  For polls = 1 To 65000
    Call SyncLog(False)
    code = ExecScalar("if test -f " & remoteDir & "/" & name & ".exit; then cat " & remoteDir & "/" & name & ".exit; else printf pending; fi")
    Call Save("heartbeat.txt", CStr(Now))
    If code <> "pending" Then
      Call SyncLog(True)
      If code <> "0" Then Call Fail("Stage " & name & " failed (" & code & "). Inspect the log. Submitted jobs are not cancelled.")
      Exit Sub
    End If
    crt.Sleep 120000
  Next
  Call Fail("Monitoring timeout; remote jobs are not cancelled.")
End Sub

Sub Main
  Call Save("state.txt", "connecting|admin2 -> fat2")
  Set tab = crt.Session.ConnectInTab("/S admin2")
  tab.Screen.Synchronous = True
  tab.Screen.IgnoreEscape = True
  Call PromptReady()
  Call QuietTerminal()
  Call Status("Connected to admin2")
  tab.Screen.Send "ssh fat2" & vbCr
  Call PromptReady()
  Call QuietTerminal()
  Call Status("Connected to fat2")
  Call ExecChecked("test ""$(id -un)"" = yzhao && test ""$(hostname -s)"" = fat2 && test -d {BASE}")

  ' Before any upload/lock/analysis, check whether the submitted phenotype
  ' filename already exists in the target project's phenotype directory.
  ' If it does, pause here and let the web page ask the user once.
  Dim duplicatePhenotype
  duplicatePhenotype = ExecScalar("if test -e " & ShellQuote({vb(BASE + "/" + folder + "/phenotype/" + filename)}) & "; then printf yes; else printf no; fi")
  If duplicatePhenotype = "yes" And Not fso.FileExists(localDir & "\\overwrite-confirmed.txt") Then
    Call Save("state.txt", "confirm|同名表型文件已存在：" & {vb(filename)} & "。请在网页确认“继续并覆盖”或“停止本次分析”。")
    tab.Screen.Send "exit" & vbCr
    Call PromptReady()
    tab.Session.Disconnect
    Exit Sub
  End If

  Call ExecChecked("mkdir -p {remote}")
  Call Status("Uploading workflow and phenotype/control files")
{uploads}
  Call Status("Starting prepare stage")
  Call Stage("prepare", {vb(launch['prepare'])})
  Call Status("Prepare stage completed")
  Call Save("state.txt", "connecting|fat2 -> 10.10.1.20")
  tab.Screen.Send "ssh 10.10.1.20" & vbCr
  Call PromptReady()
  Call QuietTerminal()
  Call Status("Connected to 10.10.1.20")
  Call ExecChecked("hostname -I | tr ' ' '\\n' | grep -Fx 10.10.1.20")
  Dim entered, archiveName
  entered = Now
  archiveName = {vb(Path(filename).stem + '_')} & Year(entered) & Right("0" & Month(entered), 2) & Right("0" & Day(entered), 2) & "-" & Right("0" & Hour(entered), 2) & Right("0" & Minute(entered), 2)
  Call Save("archive-name.txt", archiveName)
  Call ExecChecked("test -s " & remoteDir & "/archive-name.txt || printf '%s' " & ShellQuote(archiveName) & " > " & remoteDir & "/archive-name.txt")
  tab.Screen.Send {vb('cd -- ' + q(BASE + '/' + folder))} & vbCr
  Call PromptReady()
  Call ExecChecked({vb('test "$PWD" = ' + q(BASE + '/' + folder) + ' && pwd')})
  Call Status("Submitting and monitoring GWAS/plot jobs")
  Call Stage("admin", {vb(launch['admin'])})
  Call Status("GWAS and plot stages completed")
  tab.Screen.Send "exit" & vbCr
  Call PromptReady()
  Call Status("Starting final result extraction and archive")
  Call Stage("final", {vb(launch['final'])})
  Call Save("state.txt", "complete|GWAS analysis completed")
  tab.Screen.Send "printf '\\nGWAS analysis completed successfully.\\n'" & vbCr
  Call PromptReady()
  tab.Screen.Send "exit" & vbCr
  Call PromptReady()
  tab.Screen.Send "exit" & vbCr
End Sub
'''


def recovery_script(run, meta):
    """Reconnect to the same batch; never replay a submitted admin phase."""
    pop, folder, filename = meta['population'], meta['directory'], meta['filename']
    s = bridge_script(run, pop, folder, filename, meta['label'], meta.get('dataSource'))
    remote = BASE + '/.gwas-web/' + run.name
    work = BASE + '/' + folder
    batch = work + '/.gwas-runs/' + run.name
    q = shlex.quote
    begin = s.index('  Call ExecChecked("mkdir -p ', s.index('Sub Main'))
    end = s.index('  Call Save("state.txt", "connecting|fat2 ->', begin)
    # Keep the running workflow untouched; the new copy is used only for unstarted stages.
    payload = base64.b64encode(WORKFLOW.read_bytes().replace(b'\r\n', b'\n')).decode()
    path = remote + '/recovery-workflow.sh'
    commands = ['test "$(cat ' + q(remote+'/prepare.exit') + ')" = 0', ': > '+q(path+'.b64')]
    commands += ["printf '%s' '" + payload[i:i+3000] + "' >> " + q(path+'.b64') for i in range(0,len(payload),3000)]
    commands += ['base64 -d '+q(path+'.b64')+' > '+q(path)]
    s = s[:begin] + '\n'.join('  Call ExecChecked('+vb(c)+')' for c in commands) + '\n' + s[end:]
    s = s.replace(remote+'/workflow.sh', path)
    start = s.index('  Call Stage("admin",', s.index('Sub Main'))
    end = s.index('\n', start)
    original = s[start:end]
    s = s[:start] + '''  Dim recoveredCode
  recoveredCode = ExecScalar("if test -f " & remoteDir & "/admin.exit; then cat " & remoteDir & "/admin.exit; else printf pending; fi")
  If recoveredCode = "pending" Then
    If ExecScalar(''' + vb('if test -e '+q(batch+'/admin.started')+'; then printf started; else printf new; fi') + ''') = "started" Then
      Call WaitStage("admin")
    Else
''' + original + '''
    End If
  ElseIf recoveredCode <> "0" Then
    Call SyncLog(True)
    Call Fail("Remote admin phase failed with exit code " & recoveredCode & ". Review the last ERROR in the log. No jobs were resubmitted.")
  End If''' + s[end:]
    start = s.index('  Call Stage("final",', s.index('Sub Main'))
    end = s.index('\n', start)
    original = s[start:end]
    s = s[:start] + '''  recoveredCode = ExecScalar("if test -f " & remoteDir & "/final.exit; then cat " & remoteDir & "/final.exit; else printf pending; fi")
  If recoveredCode = "pending" Then
    Call ExecChecked(''' + vb('test ! -e '+q(batch+'/final.started')) + ''')
''' + original + '''
  ElseIf recoveredCode <> "0" Then
    Call SyncLog(True)
    Call Fail("Final phase needs manual repair. No analysis was replayed.")
  End If''' + s[end:]
    return s


def read_bridge(path):
    # bridge.vbs overwrites very small status/log files in place.  Windows can
    # briefly deny a simultaneous read while that handle is open, so retry the
    # read instead of letting /api/runs silently lose the run.
    last_error = None
    for _ in range(20):
        try:
            raw = path.read_bytes()
            return raw.decode('utf-16' if raw.startswith((b'\xff\xfe', b'\xfe\xff')) else 'utf-8', errors='replace')
        except OSError as exc:
            last_error = exc
            time.sleep(0.05)
    raise last_error



def _powershell_encoded(script):
    """Encode PowerShell text for -EncodedCommand (UTF-16LE)."""
    return base64.b64encode(script.encode('utf-16le')).decode('ascii')


def windows_notify(title, message, kind='info'):
    """Show a non-blocking Windows tray notification using only built-in .NET."""
    icon = 'Warning' if kind in ('warning', 'error') else 'Information'
    sound = 'Exclamation' if kind in ('warning', 'error') else 'Asterisk'

    # PowerShell single-quoted literals escape a quote by doubling it.
    ps_title = str(title).replace("'", "''").replace('\r', ' ').replace('\n', ' ')
    ps_message = str(message).replace("'", "''").replace('\r', ' ').replace('\n', ' | ')

    script = f"""
Add-Type -AssemblyName System.Windows.Forms
Add-Type -AssemblyName System.Drawing
$notify = New-Object System.Windows.Forms.NotifyIcon
$notify.Icon = [System.Drawing.SystemIcons]::{icon}
$notify.Visible = $true
$notify.BalloonTipTitle = '{ps_title}'
$notify.BalloonTipText = '{ps_message}'
$notify.BalloonTipIcon = [System.Windows.Forms.ToolTipIcon]::{icon}
[System.Media.SystemSounds]::{sound}.Play()
$notify.ShowBalloonTip(8000)
Start-Sleep -Seconds 9
$notify.Dispose()
"""
    try:
        subprocess.Popen(
            [
                'powershell.exe',
                '-NoProfile',
                '-STA',
                '-WindowStyle', 'Hidden',
                '-EncodedCommand', _powershell_encoded(script),
            ],
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0),
        )
        return True
    except OSError:
        return False


def notify_run_state(run, meta, status, detail):
    """Notify once per run and state. Marker files persist across server restarts."""
    if status not in ('complete', 'failed', 'attention'):
        return

    marker = run / ('notified-' + status + '.txt')
    with NOTIFY_LOCK:
        if marker.exists():
            return

        population = meta.get('population', 'GWAS')
        filename = meta.get('filename', '')
        archive_file = run / 'archive-name.txt'
        archive = read_bridge(archive_file).strip() if archive_file.exists() else ''

        if status == 'complete':
            title = 'GWAS 分析完成'
            parts = [f'{population} · {filename}', '全部分析流程已完成']
            if archive:
                parts.append('归档：' + archive)
            kind = 'info'
        elif status == 'failed':
            title = 'GWAS 分析失败'
            parts = [f'{population} · {filename}', detail[:180] or '请查看运行日志']
            kind = 'error'
        else:
            title = 'GWAS 需要处理'
            parts = [f'{population} · {filename}', detail[:180] or '请查看运行日志和 SecureCRT']
            kind = 'warning'

        if windows_notify(title, '\n'.join(parts), kind):
            marker.write_text(
                time.strftime('%Y-%m-%d %H:%M:%S') + '\n' + title + '\n' + '\n'.join(parts) + '\n',
                encoding='utf-8'
            )

def snapshot(run):
    meta = json.loads((run / 'meta.json').read_text('utf-8'))
    state = run / 'state.txt'
    value = read_bridge(state) if state.exists() else 'connecting|Waiting for SecureCRT'
    status, _, detail = value.partition('|')
    if meta.get('remoteDirectoryDeleted'):
        status, detail = 'attention', '服务器工作目录已由用户删除，本批次无法恢复。旧队列是否结束尚待核实；请勿复用原目录重提。新分析请填写新的工作目录，本地日志可继续保存在原文件夹，重名时自动加时间后缀。'
    if status == 'attention' and 'timeout' in detail.lower():
        detail = 'SecureCRT 通信超时，尚不能判断远程任务是否结束。关闭旧错误弹窗后，可点击“重新连接并核实 / 继续”。如果工作目录已删除，不能恢复；请待旧队列结束后重新建任务。'
        last_command = run/'last-command.txt'
        if last_command.exists():
            detail += '\n超时前操作：' + read_bridge(last_command)
    heartbeat = run / 'heartbeat.txt'
    latest = max((p.stat().st_mtime for p in (state, heartbeat) if p.exists()), default=run.stat().st_mtime)
    if status in ('connecting', 'running') and time.time() - latest > 360:
        status, detail = 'attention', '连接或监控超过6分钟未更新；远程任务可能仍在运行，请检查 SecureCRT，勿重复提交。'
    remote = run / 'remote.log'
    log = read_bridge(remote) if remote.exists() else ''
    header = meta['header']
    # Preserve Unicode terminal output; the user-facing log is UTF-8.
    try:
        Path(meta['logPath']).write_text(header + '\n' + log + '\n状态: ' + status + ' | ' + detail + '\n', 'utf-8')
    except OSError as exc:
        detail += '；本地日志写入失败: ' + str(exc)

    # mirror_logs() calls snapshot every 5 seconds even if the browser is closed,
    # so Windows notifications do not depend on the webpage being in the foreground.
    try:
        notify_run_state(run, meta, status, detail)
    except (OSError, ValueError):
        pass

    return dict(meta, status=status, detail=detail, log=log[-100000:], overwriteConfirmed=(run/'overwrite-confirmed.txt').exists())


def mirror_logs():
    while True:
        for run in RUNS.iterdir():
            if (run / 'meta.json').exists():
                try:
                    snapshot(run)
                except (OSError, ValueError):
                    pass
        time.sleep(5)


class Handler(BaseHTTPRequestHandler):
    def log_message(self, format, *args):
        # The browser polls every five seconds; keep the background log small.
        if '/api/runs' not in str(args) and '/api/health' not in str(args):
            super().log_message(format, *args)

    def send(self, data, status=200, content_type='application/json; charset=utf-8'):
        raw = json.dumps(data, ensure_ascii=False).encode() if isinstance(data, (dict, list)) else data
        self.send_response(status)
        self.send_header('Content-Type', content_type)
        self.send_header('Content-Length', str(len(raw)))
        self.send_header('Cache-Control', 'no-store')
        self.send_header('X-Content-Type-Options', 'nosniff')
        self.send_header('Content-Security-Policy', "default-src 'self'; style-src 'self' 'unsafe-inline'; script-src 'self'; connect-src 'self'; img-src 'self' data:; frame-ancestors 'none'")
        self.end_headers()
        self.wfile.write(raw)

    def send_attachment(self, raw, filename, content_type='text/plain; charset=utf-8'):
        safe_name = re.sub(r'[^A-Za-z0-9_.-]+', '_', Path(filename).name)
        self.send_response(200)
        self.send_header('Content-Type', content_type)
        self.send_header('Content-Length', str(len(raw)))
        self.send_header('Content-Disposition', f'attachment; filename="{safe_name}"')
        self.send_header('Cache-Control', 'no-store')
        self.send_header('X-Content-Type-Options', 'nosniff')
        self.end_headers()
        self.wfile.write(raw)

    def valid_host(self):
        return self.headers.get('Host') in (f'127.0.0.1:{PORT}', f'localhost:{PORT}')

    def do_GET(self):
        if not self.valid_host():
            return self.send({'error': 'Invalid host'}, 403)

        parsed = urlparse(self.path)
        if parsed.path == '/api/log-download':
            rid = parse_qs(parsed.query).get('id', [''])[0]
            if not re.fullmatch(r'[A-Za-z0-9_-]+', rid):
                return self.send({'error': '无效任务编号'}, 400)
            run = RUNS / rid
            meta_path = run / 'meta.json'
            if not meta_path.exists():
                return self.send({'error': '任务不存在'}, 404)
            try:
                # Refresh the mirrored local log before downloading it.
                snapshot(run)
                meta = json.loads(meta_path.read_text('utf-8'))
                log_path = Path(meta['logPath'])
                if not log_path.exists() or not log_path.is_file():
                    return self.send({'error': '日志文件不存在'}, 404)
                raw = log_path.read_bytes()
                return self.send_attachment(raw, log_path.name)
            except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
                return self.send({'error': '日志下载失败: ' + str(exc)}, 400)

        if self.path == '/api/health':
            return self.send({'app': 'gwas-local-console', 'version': 32})
        if self.path == '/api/config':
            return self.send({'token': TOKEN, 'securecrt': CRT.exists(), 'defaultLogDirectory': str(ROOT.parent), 'bookmarkUrl': f'http://127.0.0.1:{PORT}/', 'version': 32, 'dataSources': {p: data_source({}, p) for p in ('GPall', 'GPallInd', 'GPallJap')}})
        if self.path == '/api/runs':
            result = []
            for run in sorted(RUNS.iterdir(), reverse=True):
                meta_path = run / 'meta.json'
                if not meta_path.exists():
                    continue
                try:
                    result.append(snapshot(run))
                except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
                    # Never make a historical run disappear just because one
                    # local status/log read failed.  Return its persisted meta
                    # and make the local problem visible in the UI.
                    try:
                        meta = json.loads(meta_path.read_text('utf-8'))
                    except Exception:
                        meta = {
                            'id': run.name,
                            'population': '',
                            'directory': '',
                            'filename': '',
                            'label': '',
                            'logPath': ''
                        }
                    result.append(dict(
                        meta,
                        id=meta.get('id', run.name),
                        status='attention',
                        detail='Local snapshot error: ' + str(exc),
                        log='',
                        overwriteConfirmed=(run/'overwrite-confirmed.txt').exists()
                    ))
            return self.send(result)
        path = {'/': 'index.html', '/app.js': 'app.js', '/style.css': 'style.css'}.get(self.path)
        if not path:
            return self.send({'error': 'Not found'}, 404)
        mime = {'html': 'text/html', 'js': 'text/javascript', 'css': 'text/css'}[path.split('.')[-1]]
        raw = (ROOT / 'web' / path).read_bytes()
        if path == 'app.js':
            raw += ';\n// ---- GWAS v25 UI: robust remote-directory field + generated full path ----\n(() => {\n  const POPS = [\'GPallInd\', \'GPallJap\', \'GPall\'];\n  const BASE = \'/data9/home/yzhao/GWAS_IRGSP1.0/\';\n\n  function inferPopulation(el) {\n    if (!el) return \'\';\n    if (POPS.includes(el.value)) return el.value;\n\n    let node = el;\n    for (let n = 0; n < 7 && node; n++, node = node.parentElement) {\n      const text = (node.innerText || node.textContent || \'\').replace(/\\s+/g, \' \');\n      if (/\\bGPallInd\\b/.test(text)) return \'GPallInd\';\n      if (/\\bGPallJap\\b/.test(text)) return \'GPallJap\';\n      if (/\\bGPall\\b/.test(text)) return \'GPall\';\n    }\n    return \'\';\n  }\n\n  function selectedPopulation() {\n    const checked = [...document.querySelectorAll(\'input[type="radio"]:checked\')];\n    for (const r of checked) {\n      const pop = inferPopulation(r);\n      if (pop) return pop;\n    }\n    return \'\';\n  }\n\n  function isTextInput(x) {\n    if (!(x instanceof HTMLInputElement)) return false;\n    const t = (x.getAttribute(\'type\') || \'text\').toLowerCase();\n    return t === \'text\' || t === \'search\' || t === \'\';\n  }\n\n  function directText(el) {\n    if (!el) return \'\';\n    return [...el.childNodes]\n      .filter(n => n.nodeType === Node.TEXT_NODE)\n      .map(n => n.textContent || \'\')\n      .join(\' \')\n      .replace(/\\s+/g, \' \')\n      .trim();\n  }\n\n  function directoryInput() {\n    // 1) Find the visible label itself, then the nearest text input below/inside it.\n    const labels = [...document.querySelectorAll(\'label, h1, h2, h3, h4, h5, strong, div, p, span\')]\n      .filter(el => {\n        const d = directText(el);\n        return d.includes(\'远程目录名称\');\n      });\n\n    const textInputs = [...document.querySelectorAll(\'input\')].filter(isTextInput);\n\n    for (const label of labels) {\n      // First try a reasonably small ancestor containing one text input.\n      let node = label.parentElement;\n      for (let depth = 0; depth < 4 && node; depth++, node = node.parentElement) {\n        const found = [...node.querySelectorAll(\'input\')].filter(isTextInput);\n        if (found.length === 1) return found[0];\n      }\n\n      // Fallback: nearest text input appearing after this label in the document.\n      const after = textInputs\n        .filter(inp => label.compareDocumentPosition(inp) & Node.DOCUMENT_POSITION_FOLLOWING);\n      if (after.length) return after[0];\n    }\n\n    // 2) Existing/remembered directory values.\n    return textInputs.find(x =>\n      /^(?:(?:GPallInd|GPallJap|GPall)?miss20_)?[A-Za-z0-9][A-Za-z0-9_.-]*$/\n        .test((x.value || \'\').trim())\n      && /(?:UAV|lines|Re-fitting|miss20)/i.test(x.value || \'\')\n    ) || null;\n  }\n\n  function normalizeSuffix(value) {\n    return String(value || \'\')\n      .trim()\n      .replace(/^(?:GPallInd|GPallJap|GPall)?miss20_/, \'\');\n  }\n\n  function setNativeInputValue(input, value) {\n    if (!input || input.value === value) return;\n    const setter = Object.getOwnPropertyDescriptor(\n      window.HTMLInputElement.prototype, \'value\'\n    )?.set;\n    if (setter) setter.call(input, value);\n    else input.value = value;\n    input.dispatchEvent(new Event(\'input\', {bubbles: true}));\n    input.dispatchEvent(new Event(\'change\', {bubbles: true}));\n  }\n\n  function generatedFullPath() {\n    const input = directoryInput();\n    const pop = selectedPopulation();\n    if (!input || !pop) return \'\';\n    const suffix = normalizeSuffix(input.value);\n    if (!suffix) return BASE + pop + \'miss20_\';\n    return BASE + pop + \'miss20_\' + suffix;\n  }\n\n  function pathDisplayNode(input) {\n    if (!input) return null;\n\n    // Prefer an already existing path display near the remote-directory input.\n    const all = [...document.querySelectorAll(\'div, p, span, small\')].filter(el => {\n      if (el.children.length) return false;\n      return (el.textContent || \'\').trim().startsWith(BASE);\n    });\n    if (!all.length) return null;\n\n    const ir = input.getBoundingClientRect();\n    all.sort((a, b) => {\n      const ar = a.getBoundingClientRect();\n      const br = b.getBoundingClientRect();\n      const da = Math.abs(ar.top - ir.bottom) + Math.abs(ar.left - ir.left) * 0.1;\n      const db = Math.abs(br.top - ir.bottom) + Math.abs(br.left - ir.left) * 0.1;\n      return da - db;\n    });\n    return all[0];\n  }\n\n  function refreshDirectory() {\n    const input = directoryInput();\n    if (!input) return;\n\n    // v25 field contains ONLY what comes after "miss20_".\n    const suffix = normalizeSuffix(input.value);\n    if (suffix !== input.value) {\n      setNativeInputValue(input, suffix);\n    }\n\n    const full = generatedFullPath();\n    const display = pathDisplayNode(input);\n    if (display && full) display.textContent = full;\n  }\n\n  // Population change: only the generated full path changes.\n  document.addEventListener(\'change\', ev => {\n    if (ev.target instanceof HTMLInputElement && ev.target.type === \'radio\') {\n      if (inferPopulation(ev.target)) setTimeout(refreshDirectory, 0);\n    }\n  }, true);\n\n  document.addEventListener(\'click\', ev => {\n    const card = ev.target.closest(\'label, [role="radio"]\');\n    if (!card) return;\n    const text = (card.innerText || card.textContent || \'\').replace(/\\s+/g, \' \');\n    if (/\\b(?:GPallInd|GPallJap|GPall)\\b/.test(text)) {\n      setTimeout(refreshDirectory, 0);\n    }\n  }, true);\n\n  // Typing the reusable suffix updates the full path immediately.\n  document.addEventListener(\'input\', ev => {\n    const input = directoryInput();\n    if (input && ev.target === input) {\n      const suffix = normalizeSuffix(input.value);\n      if (suffix !== input.value) setNativeInputValue(input, suffix);\n      setTimeout(refreshDirectory, 0);\n    }\n  }, true);\n\n  // Preserve the run selector population prefix enhancement.\n  let updatingRuns = false;\n  async function prefixRunLabels() {\n    if (updatingRuns) return;\n    updatingRuns = true;\n    try {\n      const runs = await fetch(\'/api/runs\', {cache: \'no-store\'})\n        .then(r => r.ok ? r.json() : []);\n      const byId = new Map(runs.map(r => [String(r.id), r]));\n      for (const select of document.querySelectorAll(\'select\')) {\n        for (const option of select.options) {\n          const run = byId.get(String(option.value));\n          if (!run || !run.population) continue;\n          const prefix = run.population + \'_\';\n          if (!option.textContent.startsWith(prefix)) {\n            option.textContent = prefix + option.textContent;\n          }\n        }\n      }\n    } catch (_) {\n    } finally {\n      updatingRuns = false;\n    }\n  }\n\n  function refreshAll() {\n    refreshDirectory();\n    prefixRunLabels();\n  }\n\n  window.addEventListener(\'DOMContentLoaded\', () => setTimeout(refreshAll, 80));\n  setTimeout(refreshAll, 350);\n\n  const observer = new MutationObserver(() => setTimeout(refreshAll, 0));\n  observer.observe(document.body, {childList: true, subtree: true});\n})();'.encode('utf-8')
        if path == 'app.js':
            raw += ";\n// ---- GWAS v28: download the full local run log with the real log filename ----\n(() => {\n  async function downloadSelectedRunLog(ev) {\n    const target = ev.target.closest('button, a');\n    if (!target) return;\n    const label = (target.innerText || target.textContent || '').replace(/\\s+/g, ' ').trim();\n    if (!label.includes('下载日志')) return;\n\n    const runs = await fetch('/api/runs', {cache: 'no-store'})\n      .then(r => r.ok ? r.json() : []);\n\n    const byId = new Map(runs.map(r => [String(r.id), r]));\n    let runId = '';\n\n    // The run selector uses the run ID as option value.\n    for (const select of document.querySelectorAll('select')) {\n      if (byId.has(String(select.value))) {\n        runId = String(select.value);\n        break;\n      }\n    }\n\n    if (!runId) return; // fall back to the original page behavior if no run can be resolved\n\n    ev.preventDefault();\n    ev.stopImmediatePropagation();\n\n    // The server sends Content-Disposition using the actual log basename,\n    // e.g. SJ2024_GWASrun_GPall.log.\n    window.location.href = '/api/log-download?id=' + encodeURIComponent(runId);\n  }\n\n  document.addEventListener('click', ev => {\n    // Keep this handler synchronous enough to intercept the original click.\n    const target = ev.target.closest('button, a');\n    if (!target) return;\n    const label = (target.innerText || target.textContent || '').replace(/\\s+/g, ' ').trim();\n    if (!label.includes('下载日志')) return;\n\n    ev.preventDefault();\n    ev.stopImmediatePropagation();\n\n    fetch('/api/runs', {cache: 'no-store'})\n      .then(r => r.ok ? r.json() : [])\n      .then(runs => {\n        const byId = new Map(runs.map(r => [String(r.id), r]));\n        let runId = '';\n        for (const select of document.querySelectorAll('select')) {\n          if (byId.has(String(select.value))) {\n            runId = String(select.value);\n            break;\n          }\n        }\n        if (runId) {\n          window.location.href = '/api/log-download?id=' + encodeURIComponent(runId);\n        }\n      })\n      .catch(() => {});\n  }, true);\n})();".encode('utf-8')
        if path == 'app.js':
            raw += ';\n// ---- GWAS v29: one-time confirmation when remote phenotype filename already exists ----\n(() => {\n  let modalRunId = \'\';\n  let checking = false;\n\n  async function apiPost(path, payload) {\n    const cfg = await fetch(\'/api/config\', {cache: \'no-store\'}).then(r => r.json());\n    const r = await fetch(path, {\n      method: \'POST\',\n      headers: {\'Content-Type\': \'application/json\', \'X-GWAS-Token\': cfg.token},\n      body: JSON.stringify(payload)\n    });\n    const body = await r.json().catch(() => ({}));\n    if (!r.ok) throw new Error(body.error || (\'HTTP \' + r.status));\n    return body;\n  }\n\n  function removeModal() {\n    const old = document.getElementById(\'gwas-overwrite-confirm-modal\');\n    if (old) old.remove();\n    modalRunId = \'\';\n  }\n\n  function showModal(run) {\n    if (!run || !run.id || modalRunId === String(run.id)) return;\n    removeModal();\n    modalRunId = String(run.id);\n\n    const overlay = document.createElement(\'div\');\n    overlay.id = \'gwas-overwrite-confirm-modal\';\n    Object.assign(overlay.style, {\n      position: \'fixed\', inset: \'0\', zIndex: \'2147483647\',\n      background: \'rgba(0,0,0,.42)\', display: \'flex\',\n      alignItems: \'center\', justifyContent: \'center\', padding: \'24px\'\n    });\n\n    const box = document.createElement(\'div\');\n    Object.assign(box.style, {\n      width: \'min(620px, 92vw)\', background: \'#fff\', borderRadius: \'12px\',\n      boxShadow: \'0 18px 60px rgba(0,0,0,.28)\', padding: \'24px\',\n      fontFamily: \'Arial, sans-serif\', color: \'#111\'\n    });\n\n    const title = document.createElement(\'div\');\n    title.textContent = \'发现同名表型文件\';\n    Object.assign(title.style, {fontSize: \'20px\', fontWeight: \'700\', marginBottom: \'14px\'});\n\n    const text = document.createElement(\'div\');\n    text.innerHTML =\n      \'远程 <b>phenotype</b> 目录中已经存在同名文件：<br>\' +\n      \'<code style="display:inline-block;margin:8px 0;padding:4px 7px;background:#f5f5f5;border-radius:5px;">\' +\n      String(run.filename || \'\') +\n      \'</code><br><br>\' +\n      \'<br><br><b>请确认本次表型文件名称和内容是否正确。</b><br>\' +\n      \'选择“继续并覆盖”后，本次运行会覆盖之前的结果。\';\n\n    const buttons = document.createElement(\'div\');\n    Object.assign(buttons.style, {\n      display: \'flex\', justifyContent: \'flex-end\', gap: \'12px\', marginTop: \'22px\'\n    });\n\n    const stop = document.createElement(\'button\');\n    stop.type = \'button\';\n    stop.textContent = \'停止本次分析\';\n    Object.assign(stop.style, {\n      padding: \'9px 16px\', border: \'1px solid #bbb\', borderRadius: \'7px\',\n      background: \'#fff\', color: \'#111827\', cursor: \'pointer\', fontWeight: \'600\'\n    });\n\n    const go = document.createElement(\'button\');\n    go.type = \'button\';\n    go.textContent = \'继续并覆盖\';\n    Object.assign(go.style, {\n      padding: \'9px 16px\', border: \'0\', borderRadius: \'7px\',\n      background: \'#111827\', color: \'#fff\', cursor: \'pointer\', fontWeight: \'600\'\n    });\n\n    async function decide(proceed) {\n      stop.disabled = true;\n      go.disabled = true;\n      try {\n        await apiPost(\'/api/overwrite-decision\', {id: run.id, proceed});\n        removeModal();\n      } catch (e) {\n        stop.disabled = false;\n        go.disabled = false;\n        alert(\'提交确认失败：\' + e.message);\n      }\n    }\n\n    stop.addEventListener(\'click\', () => decide(false));\n    go.addEventListener(\'click\', () => decide(true));\n\n    buttons.append(stop, go);\n    box.append(title, text, buttons);\n    overlay.append(box);\n    document.body.append(overlay);\n  }\n\n  async function checkOverwriteConfirmation() {\n    if (checking) return;\n    checking = true;\n    try {\n      const runs = await fetch(\'/api/runs\', {cache: \'no-store\'}).then(r => r.ok ? r.json() : []);\n      const pending = runs.find(r => r && r.status === \'confirm\');\n      if (pending) {\n        showModal(pending);\n      } else if (modalRunId) {\n        removeModal();\n      }\n    } catch (_) {\n    } finally {\n      checking = false;\n    }\n  }\n\n  setInterval(checkOverwriteConfirmation, 1500);\n  setTimeout(checkOverwriteConfirmation, 300);\n})();'.encode('utf-8')
        if path == 'app.js':
            raw += ';\n// ---- GWAS v31: keep the Data Version field aligned with the selected population ----\n(() => {\n  const POPS = [\'GPallInd\', \'GPallJap\', \'GPall\'];\n  const AUTO_RE = /^(?:miss20-v1|(?:GPallInd|GPallJap|GPall)miss20-864lines)$/;\n\n  function inferPopulation(el) {\n    if (!el) return \'\';\n    if (POPS.includes(el.value)) return el.value;\n    let node = el;\n    for (let n = 0; n < 7 && node; n++, node = node.parentElement) {\n      const t = (node.innerText || node.textContent || \'\').replace(/\\s+/g, \' \');\n      if (/\\bGPallInd\\b/.test(t)) return \'GPallInd\';\n      if (/\\bGPallJap\\b/.test(t)) return \'GPallJap\';\n      if (/\\bGPall\\b/.test(t)) return \'GPall\';\n    }\n    return \'\';\n  }\n\n  function selectedPopulation() {\n    for (const r of document.querySelectorAll(\'input[type="radio"]:checked\')) {\n      const p = inferPopulation(r);\n      if (p) return p;\n    }\n    return \'\';\n  }\n\n  function directText(el) {\n    if (!el) return \'\';\n    return [...el.childNodes]\n      .filter(n => n.nodeType === Node.TEXT_NODE)\n      .map(n => n.textContent || \'\')\n      .join(\' \')\n      .replace(/\\s+/g, \' \')\n      .trim();\n  }\n\n  function versionInput() {\n    const labels = [...document.querySelectorAll(\'label, strong, div, p, span\')]\n      .filter(el => directText(el).includes(\'数据版本\'));\n    for (const label of labels) {\n      let node = label.parentElement;\n      for (let depth = 0; depth < 4 && node; depth++, node = node.parentElement) {\n        const inputs = [...node.querySelectorAll(\'input\')].filter(x => {\n          const type = (x.getAttribute(\'type\') || \'text\').toLowerCase();\n          return type === \'text\' || type === \'\';\n        });\n        if (inputs.length === 1) return inputs[0];\n      }\n    }\n    return null;\n  }\n\n  function setNativeValue(input, value) {\n    if (!input || input.value === value) return;\n    const setter = Object.getOwnPropertyDescriptor(\n      window.HTMLInputElement.prototype, \'value\'\n    )?.set;\n    if (setter) setter.call(input, value);\n    else input.value = value;\n    input.dispatchEvent(new Event(\'input\', {bubbles: true}));\n    input.dispatchEvent(new Event(\'change\', {bubbles: true}));\n  }\n\n  function refreshVersion(forceAuto = false) {\n    const input = versionInput();\n    const pop = selectedPopulation();\n    if (!input || !pop) return;\n    const current = String(input.value || \'\').trim();\n    const desired = pop + \'miss20-864lines\';\n\n    // Blank values and older automatic defaults are repaired automatically.\n    // A genuinely custom version entered by the user is left untouched.\n    if (!current || AUTO_RE.test(current) || forceAuto) {\n      setNativeValue(input, desired);\n    }\n  }\n\n  document.addEventListener(\'change\', ev => {\n    if (ev.target instanceof HTMLInputElement && ev.target.type === \'radio\') {\n      if (inferPopulation(ev.target)) setTimeout(() => refreshVersion(false), 0);\n    }\n  }, true);\n\n  document.addEventListener(\'click\', ev => {\n    const card = ev.target.closest(\'label, [role="radio"]\');\n    if (!card) return;\n    const t = (card.innerText || card.textContent || \'\').replace(/\\s+/g, \' \');\n    if (/\\b(?:GPallInd|GPallJap|GPall)\\b/.test(t)) {\n      setTimeout(() => refreshVersion(false), 0);\n    }\n  }, true);\n\n  setTimeout(() => refreshVersion(false), 100);\n  setTimeout(() => refreshVersion(false), 500);\n\n  const observer = new MutationObserver(() => setTimeout(() => refreshVersion(false), 0));\n  observer.observe(document.body, {childList: true, subtree: true});\n})();'.encode('utf-8')
        if path == 'app.js':
            raw += ';\n// ---- GWAS v32: robust data-source state + visible submit errors ----\n(() => {\n  const POPS = [\'GPallInd\', \'GPallJap\', \'GPall\'];\n  const AUTO_RE = /^(?:miss20-v1|(?:GPallInd|GPallJap|GPall)miss20-864lines)$/;\n  let defaults = null;\n  let toastTimer = null;\n\n  function showToast(message, kind = \'error\') {\n    let box = document.getElementById(\'gwas-v32-toast\');\n    if (!box) {\n      box = document.createElement(\'div\');\n      box.id = \'gwas-v32-toast\';\n      Object.assign(box.style, {\n        position: \'fixed\', right: \'22px\', top: \'22px\', zIndex: \'2147483647\',\n        maxWidth: \'620px\', padding: \'13px 16px\', borderRadius: \'9px\',\n        boxShadow: \'0 10px 35px rgba(0,0,0,.22)\',\n        fontFamily: \'Arial, sans-serif\', fontSize: \'14px\',\n        whiteSpace: \'pre-wrap\', lineHeight: \'1.45\'\n      });\n      document.body.appendChild(box);\n    }\n    box.style.background = kind === \'ok\' ? \'#ecfdf5\' : \'#fff1f2\';\n    box.style.color = kind === \'ok\' ? \'#065f46\' : \'#991b1b\';\n    box.style.border = kind === \'ok\' ? \'1px solid #a7f3d0\' : \'1px solid #fecdd3\';\n    box.textContent = message;\n    box.style.display = \'block\';\n    clearTimeout(toastTimer);\n    toastTimer = setTimeout(() => { box.style.display = \'none\'; }, kind === \'ok\' ? 3500 : 10000);\n  }\n\n  function inferPopulation(el) {\n    if (!el) return \'\';\n    if (POPS.includes(el.value)) return el.value;\n    let node = el;\n    for (let n = 0; n < 7 && node; n++, node = node.parentElement) {\n      const t = (node.innerText || node.textContent || \'\').replace(/\\s+/g, \' \');\n      if (/\\bGPallInd\\b/.test(t)) return \'GPallInd\';\n      if (/\\bGPallJap\\b/.test(t)) return \'GPallJap\';\n      if (/\\bGPall\\b/.test(t)) return \'GPall\';\n    }\n    return \'\';\n  }\n\n  function selectedPopulation() {\n    for (const r of document.querySelectorAll(\'input[type="radio"]:checked\')) {\n      const p = inferPopulation(r);\n      if (p) return p;\n    }\n    return \'\';\n  }\n\n  function isTextInput(el) {\n    if (!(el instanceof HTMLInputElement)) return false;\n    const t = (el.getAttribute(\'type\') || \'text\').toLowerCase();\n    return t === \'text\' || t === \'search\' || t === \'\';\n  }\n\n  function visibleTextInputs() {\n    return [...document.querySelectorAll(\'input\')].filter(el => {\n      if (!isTextInput(el)) return false;\n      const r = el.getBoundingClientRect();\n      return r.width > 0 && r.height > 0;\n    });\n  }\n\n  function directText(el) {\n    if (!el) return \'\';\n    return [...el.childNodes]\n      .filter(n => n.nodeType === Node.TEXT_NODE)\n      .map(n => n.textContent || \'\')\n      .join(\' \').replace(/\\s+/g, \' \').trim();\n  }\n\n  function versionInput() {\n    const inputs = visibleTextInputs();\n\n    // Strongest fallback for the current page: the Data Version input is the\n    // text input immediately before the Sample List path input.\n    const sampleIdx = inputs.findIndex(x => /\\.sampleList\\s*$/.test((x.value || \'\').trim()));\n    if (sampleIdx > 0) return inputs[sampleIdx - 1];\n\n    // Locate the exact "数据版本" label and choose the nearest following text input.\n    const labels = [...document.querySelectorAll(\'label, strong, div, p, span\')]\n      .filter(el => directText(el).includes(\'数据版本\'));\n    for (const label of labels) {\n      const lr = label.getBoundingClientRect();\n      const candidates = inputs\n        .filter(inp => label.compareDocumentPosition(inp) & Node.DOCUMENT_POSITION_FOLLOWING)\n        .map(inp => {\n          const r = inp.getBoundingClientRect();\n          return {inp, score: Math.max(0, r.top - lr.bottom) * 10 + Math.abs(r.left - lr.left)};\n        })\n        .sort((a, b) => a.score - b.score);\n      if (candidates.length) return candidates[0].inp;\n    }\n\n    // Section fallback: first text input inside the "基因型数据源" block.\n    const headings = [...document.querySelectorAll(\'div, p, span, strong, h1, h2, h3, h4\')]\n      .filter(el => directText(el).includes(\'基因型数据源\'));\n    for (const h of headings) {\n      let node = h.parentElement;\n      for (let depth = 0; depth < 6 && node; depth++, node = node.parentElement) {\n        const found = [...node.querySelectorAll(\'input\')].filter(isTextInput);\n        if (found.length >= 5) return found[0];\n      }\n    }\n    return null;\n  }\n\n  function setNativeValue(input, value) {\n    if (!input || input.value === value) return;\n    const setter = Object.getOwnPropertyDescriptor(\n      window.HTMLInputElement.prototype, \'value\'\n    )?.set;\n    if (setter) setter.call(input, value);\n    else input.value = value;\n    input.dispatchEvent(new Event(\'input\', {bubbles: true}));\n    input.dispatchEvent(new Event(\'change\', {bubbles: true}));\n  }\n\n  async function loadDefaults() {\n    if (defaults) return defaults;\n    try {\n      const cfg = await window.__gwasNativeFetch(\'/api/config\', {cache: \'no-store\'}).then(r => r.json());\n      defaults = cfg.dataSources || {};\n    } catch (_) {\n      defaults = {};\n    }\n    return defaults;\n  }\n\n  async function syncVersionField() {\n    const pop = selectedPopulation();\n    if (!pop) return;\n    const ds = await loadDefaults();\n    const desired = ds?.[pop]?.version || (pop + \'miss20-864lines\');\n    const input = versionInput();\n    if (!input) return;\n    const current = String(input.value || \'\').trim();\n    if (!current || AUTO_RE.test(current)) {\n      setNativeValue(input, desired);\n    }\n  }\n\n  // Patch the actual request body as a final guard.  This updates submitted\n  // state even if an old browser-saved empty value re-renders the field.\n  if (!window.__gwasNativeFetch) {\n    window.__gwasNativeFetch = window.fetch.bind(window);\n    window.fetch = async function(resource, init = {}) {\n      const url = typeof resource === \'string\' ? resource : (resource?.url || \'\');\n      let patchedInit = init;\n\n      if ((url.endsWith(\'/api/run\') || url.endsWith(\'/api/preview\')) &&\n          init && typeof init.body === \'string\') {\n        try {\n          const body = JSON.parse(init.body);\n          const pop = body.population;\n          if (POPS.includes(pop)) {\n            const ds = await loadDefaults();\n            const desired = ds?.[pop]?.version || (pop + \'miss20-864lines\');\n            body.dataSource = (body.dataSource && typeof body.dataSource === \'object\')\n              ? body.dataSource : {};\n            const current = String(body.dataSource.version || \'\').trim();\n            if (!current || AUTO_RE.test(current)) body.dataSource.version = desired;\n            patchedInit = {...init, body: JSON.stringify(body)};\n          }\n        } catch (_) {}\n      }\n\n      let response;\n      try {\n        response = await window.__gwasNativeFetch(resource, patchedInit);\n      } catch (err) {\n        if (url.endsWith(\'/api/run\')) {\n          showToast(\'提交失败：无法连接本地 GWAS 服务\\n\' + String(err), \'error\');\n        }\n        throw err;\n      }\n\n      if (url.endsWith(\'/api/run\')) {\n        if (!response.ok) {\n          try {\n            const body = await response.clone().json();\n            showToast(\'提交失败：\' + (body.error || (\'HTTP \' + response.status)), \'error\');\n          } catch (_) {\n            showToast(\'提交失败：HTTP \' + response.status, \'error\');\n          }\n        } else {\n          try {\n            const body = await response.clone().json();\n            showToast(\'已创建新任务：\' + (body.id || \'提交成功\'), \'ok\');\n          } catch (_) {\n            showToast(\'新任务已提交\', \'ok\');\n          }\n        }\n      }\n      return response;\n    };\n  }\n\n  document.addEventListener(\'change\', ev => {\n    if (ev.target instanceof HTMLInputElement && ev.target.type === \'radio\' && inferPopulation(ev.target)) {\n      setTimeout(syncVersionField, 0);\n    }\n  }, true);\n\n  document.addEventListener(\'click\', ev => {\n    const card = ev.target.closest(\'label, [role="radio"]\');\n    if (!card) return;\n    const t = (card.innerText || card.textContent || \'\').replace(/\\s+/g, \' \');\n    if (/\\b(?:GPallInd|GPallJap|GPall)\\b/.test(t)) setTimeout(syncVersionField, 0);\n  }, true);\n\n  setTimeout(syncVersionField, 120);\n  setTimeout(syncVersionField, 600);\n  setInterval(syncVersionField, 2500);\n})();'.encode('utf-8')
        self.send(raw, content_type=mime + '; charset=utf-8')

    def do_POST(self):
        if not self.valid_host() or self.headers.get('X-GWAS-Token') != TOKEN:
            return self.send({'error': '请求校验失败，请刷新页面'}, 403)
        if self.path not in ('/api/preview', '/api/run', '/api/recover', '/api/overwrite-decision'):
            return self.send({'error': 'Not found'}, 404)
        try:
            length = int(self.headers.get('Content-Length', '0'))
            if not 0 < length <= 29 * 1024 * 1024:
                raise ValueError('请求过大或为空')
            data = json.loads(self.rfile.read(length))
            if self.path == '/api/overwrite-decision':
                rid = data.get('id', '')
                proceed = data.get('proceed')
                if not re.fullmatch(r'[A-Za-z0-9_-]+', rid):
                    raise ValueError('无效任务编号')
                if proceed not in (True, False):
                    raise ValueError('覆盖确认参数无效')
                with LOCK:
                    run = RUNS / rid
                    meta_path = run / 'meta.json'
                    if not meta_path.exists():
                        raise ValueError('任务不存在')
                    current = snapshot(run)
                    if current['status'] != 'confirm':
                        raise ValueError('此任务当前不需要同名表型确认')
                    if not proceed:
                        (run/'state.txt').write_text(
                            'cancelled|用户选择停止：远程 phenotype 中已有同名表型文件，本次分析未进入 prepare/qsub。',
                            'utf-8'
                        )
                        return self.send(snapshot(run))
                    if not CRT.exists():
                        raise ValueError('未找到 SecureCRT: ' + str(CRT))
                    # Explicit one-run authorization. The bridge checks this marker
                    # only after it verifies the same remote phenotype collision.
                    (run/'overwrite-confirmed.txt').write_text(
                        time.strftime('%Y-%m-%d %H:%M:%S') + '\n', 'utf-8'
                    )
                    meta = json.loads(meta_path.read_text('utf-8'))
                    meta['overwriteConfirmed'] = True
                    meta_path.write_text(json.dumps(meta, ensure_ascii=False), 'utf-8')
                    (run/'state.txt').write_text(
                        'connecting|已确认同名表型；重新连接并继续，本次同名目标允许覆盖',
                        'utf-8'
                    )
                    startup = subprocess.STARTUPINFO()
                    startup.dwFlags |= subprocess.STARTF_USESHOWWINDOW
                    startup.wShowWindow = 0
                    try:
                        subprocess.Popen([str(CRT), '/SCRIPT', str(run/'bridge.vbs')], startupinfo=startup)
                    except OSError as exc:
                        (run/'state.txt').write_text('failed|SecureCRT overwrite-resume launch failed', 'utf-8')
                        raise ValueError('SecureCRT 重新启动失败: ' + str(exc)) from exc
                    return self.send(snapshot(run))
            if self.path == '/api/recover':
                rid = data.get('id', '')
                if not re.fullmatch(r'[A-Za-z0-9_-]+', rid):
                    raise ValueError('无效任务编号')
                with LOCK:
                    run = RUNS / rid
                    meta = json.loads((run/'meta.json').read_text('utf-8'))
                    if meta.get('remoteDirectoryDeleted'):
                        raise ValueError('工作目录已删除，不能恢复此批次。请新建分析配置。')
                    for prior in RUNS.iterdir():
                        if prior != run and (prior/'meta.json').exists() and snapshot(prior)['status'] in ('running', 'connecting'):
                            raise ValueError('另一个任务正在运行，请勿同时启动恢复脚本')
                    if snapshot(run)['status'] not in ('attention', 'failed'):
                        raise ValueError('仅可恢复需要处理或连接失败的任务')
                    target = run/'reconnect.vbs'
                    target.write_text(recovery_script(run, meta), 'utf-16')
                    (run/'state.txt').write_text('connecting|Reconnecting; checking remote stages without duplicate submission', 'utf-8')
                    startup = subprocess.STARTUPINFO()
                    startup.dwFlags |= subprocess.STARTF_USESHOWWINDOW
                    startup.wShowWindow = 0
                    try:
                        subprocess.Popen([str(CRT), '/SCRIPT', str(target)], startupinfo=startup)
                    except OSError:
                        (run/'state.txt').write_text('failed|SecureCRT reconnect launch failed', 'utf-8')
                        raise
                    return self.send(snapshot(run))
            pop, folder, filename, content, log = validate(data)
            label = batch_label(data, filename)
            source = data_source(data, pop)
            if self.path == '/api/preview':
                log = available_log_path(log)
                return self.send({'remotePath': BASE + '/' + folder, 'logPath': str(log), 'population': pop, 'label': label,
                                  'traits': phenotype_traits(content, pop, label),
                                  'source': source,
                                  'workflow': WORKFLOW.read_text('utf-8')})
            with LOCK:
                if not CRT.exists():
                    raise ValueError('未找到 SecureCRT: ' + str(CRT))
                for prior in RUNS.iterdir():
                    if (prior / 'meta.json').exists():
                        previous = snapshot(prior)
                        if previous['status'] in ('running', 'connecting'):
                            raise ValueError('已有任务正在运行或连接中，请等待当前任务结束后再提交。')
                log = available_log_path(log)
                rid = time.strftime('%Y%m%d-%H%M%S-') + secrets.token_hex(4)
                run = RUNS / rid
                run.mkdir()
                header = '\n'.join(['GWAS run ' + rid, '开始时间: ' + time.strftime('%Y-%m-%d %H:%M:%S'), '群体: ' + pop,
                                     '表型: ' + filename, '远程目录: ' + BASE + '/' + folder,
                                     '本次标签: ' + label, '批次标识: ' + Path(filename).stem, '仅处理本批生成清单；染色体01–12；成功后按批次归档',
                                     '数据源: ' + json.dumps(source, ensure_ascii=False), 'qsub 提交间隔: 1秒', '同名表型策略: 远程 phenotype 同名时网页确认一次；确认后仅覆盖同名目标，不做批量清理',
                                     '节点: admin2 → fat2 → 10.10.1.20 → fat2', '工作流源码:', WORKFLOW.read_text('utf-8')])
                # Exclusive create protects a log another process may have created.
                with log.open('x', encoding='utf-8') as f:
                    f.write(header)
                (run / 'phenotype.upload').write_bytes(content)
                meta = dict(id=rid, population=pop, directory=folder, filename=filename, label=label, dataSource=source, logPath=str(log), header=header)
                (run / 'meta.json').write_text(json.dumps(meta, ensure_ascii=False), 'utf-8')
                (run / 'bridge.vbs').write_text(bridge_script(run, pop, folder, filename, label, source), 'utf-16')
                (run / 'state.txt').write_text('connecting|Starting SecureCRT', 'ascii')
                try:
                    startup = subprocess.STARTUPINFO()
                    startup.dwFlags |= subprocess.STARTF_USESHOWWINDOW
                    startup.wShowWindow = 0
                    subprocess.Popen([str(CRT), '/SCRIPT', str(run / 'bridge.vbs')], startupinfo=startup)
                except Exception as exc:
                    (run / 'state.txt').write_text('failed|Cannot launch SecureCRT: ' + str(exc), 'ascii', errors='replace')
                    raise ValueError('SecureCRT 启动失败: ' + str(exc)) from exc
                return self.send(snapshot(run), 201)
        except (ValueError, OSError, TypeError) as exc:
            self.send({'error': str(exc)}, 400)


if __name__ == '__main__':
    RUNS.mkdir(exist_ok=True)
    threading.Thread(target=mirror_logs, daemon=True).start()
    print(f'GWAS console: http://127.0.0.1:{PORT}', flush=True)
    ThreadingHTTPServer(('127.0.0.1', PORT), Handler).serve_forever()
