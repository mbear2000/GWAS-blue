import json
import server
run=server.RUNS/'20260920-084444-4a8c658a'
m=json.loads((run/'meta.json').read_text('utf-8'))
s=server.bridge_script(run,m['population'],m['directory'],m['filename'],m['label'],m['dataSource'])
probe=run/'inspection'
probe.mkdir(exist_ok=True)
s=s.replace('localDir = '+server.vb(run),'localDir = '+server.vb(probe))
remote=server.BASE+'/.gwas-web/'+run.name
batch=server.BASE+'/'+m['directory']+'/.gwas-runs/'+run.name
command=f'for x in prepare admin final; do printf "%s=" "$x"; cat {remote}/$x.exit 2>/dev/null || echo pending; done; wc -l {batch}/emmax.jobs {batch}/plot.jobs; tail -c 2500 {remote}/run.log; printf "\\nQUEUE\\n"; qstat -u yzhao | tail -15'
s=s[:s.index('Sub Main')]+'''Sub Main
  Set tab = crt.Session.ConnectInTab("/S admin2")
  tab.Screen.Synchronous = True
  tab.Screen.IgnoreEscape = True
  Call PromptReady()
  tab.Screen.Send "ssh fat2" & vbCr
  Call PromptReady()
  tab.Screen.Send "ssh 10.10.1.20" & vbCr
  Call PromptReady()
  Call Save("report.txt", ExecRemote('''+server.vb(command)+'''))
  Call Save("state.txt", "complete|inspection")
End Sub
'''
(probe/'inspect.vbs').write_text(s,'utf-16')
print(probe/'inspect.vbs')
