"""Build a guarded recovery bridge for the explicitly approved prepared batch."""
import json
import subprocess
import server

run = server.RUNS / '20260920-084444-4a8c658a'
meta = json.loads((run / 'meta.json').read_text('utf-8'))
assert meta['population'] == 'GPallJap'
remote = server.BASE + '/.gwas-web/' + run.name
work = server.BASE + '/' + meta['directory']
batch = work + '/.gwas-runs/' + run.name
script = server.bridge_script(run, meta['population'], meta['directory'], meta['filename'], meta['label'], meta['dataSource'])
start = script.index('  Call ExecChecked("mkdir -p ', script.index('Sub Main'))
end = script.index('  Call Save("state.txt", "connecting|fat2 -> 10.10.1.20")', start)
q = server.shlex.quote
guard = ('test "$(cat ' + q(remote+'/prepare.exit') + ')" = 0'
         + ' && test "$(cat ' + q(work+'/.gwas-active/owner') + ')" = ' + q(run.name)
         + ' && test -s ' + q(batch+'/manifest.tsv')
         + ' && test ! -e ' + q(batch+'/admin.started')
         + ' && test ! -e ' + q(batch+'/emmax.jobs')
         + ' && test ! -e ' + q(remote+'/admin.exit')
         + ' && grep -q "sleep 5" ' + q(remote+'/workflow.sh'))
script = script[:start] + '  Call ExecChecked(' + server.vb(guard) + ')\n' + script[end:]
# Inspect queue tools and prepared inputs on the actual submission node first.
point = script.index('  Call Stage("admin",', script.index('Sub Main'))
script = script[:point] + '  Call ExecChecked("command -v qstat && command -v qsub && qstat -u yzhao")\n' + script[point:]
target = run/'resume-admin.vbs'
target.write_text(script, 'utf-16')
syntax = run/'resume-syntax-check.vbs'
syntax.write_text('\n'.join(script.replace('Option Explicit', 'Option Explicit\nDim crt').splitlines()[2:]), 'utf-16')
result = subprocess.run(['cscript.exe', '//nologo', str(syntax)], capture_output=True)
if result.returncode:
    raise RuntimeError(result.stdout + result.stderr)
print(target)
