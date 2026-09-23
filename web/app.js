const $ = id => document.getElementById(id);
let config, runs = [], selected = '', lastStatus = new Map();
const sourceKeys=['version','sample','tped','kinship','maf'];
let savedSources={};
try { savedSources=JSON.parse(localStorage.getItem('gwas-data-sources')||'{}'); } catch {}
const population=()=>document.querySelector('input[name=population]:checked').value;
const sourceValue=()=>Object.fromEntries(sourceKeys.map(k=>[k,$('source-'+k).value.trim()]));
function showSource(){const s=savedSources[population()]||config.dataSources[population()];for(const k of sourceKeys)$('source-'+k).value=s[k]||'';}
for(const k of sourceKeys)$('source-'+k).addEventListener('change',()=>{savedSources[population()]=sourceValue();try{localStorage.setItem('gwas-data-sources',JSON.stringify(savedSources))}catch{}});
for(const radio of document.querySelectorAll('input[name=population]'))radio.addEventListener('change',()=>{if(config?.dataSources)showSource()});
$('resetSource').onclick=()=>{delete savedSources[population()];try{localStorage.setItem('gwas-data-sources',JSON.stringify(savedSources))}catch{}showSource()};
const labels = {connecting:'连接中', running:'运行中', complete:'分析完成', attention:'需要处理', failed:'启动失败'};
const renderBase=render;
render=function(){renderBase();const r=runs.find(r=>r.id===selected);$('recoveryHelp').hidden=!r||!['attention','failed'].includes(r.status);$('recover').hidden=!!r?.remoteDirectoryDeleted};
$('recover').onclick=async()=>{const r=runs.find(r=>r.id===selected);if(!r)return;$('recover').disabled=true;try{await api('/api/recover',{id:r.id});message('正在重新连接，核实远程阶段；不会重复提交已启动的 GWAS。');await refresh()}catch(e){message(e.message,true)}finally{$('recover').disabled=false}};
async function api(url, data) {
  const response = await fetch(url, data ? {method:'POST',headers:{'Content-Type':'application/json','X-GWAS-Token':config.token},body:JSON.stringify(data)} : {});
  const result = await response.json();
  if (!response.ok) throw Error(result.error || '请求失败');
  return result;
}
function message(text, error=false) { $('message').textContent=text; $('message').className=error?'error':''; }
async function payload() {
  if (!$('form').reportValidity()) throw Error('请补全分析配置');
  const file=$('file').files[0];
  if(file.size > 20*1024*1024) throw Error('文件不能超过20 MB');
  const content=await new Promise((resolve,reject)=>{const r=new FileReader();r.onload=()=>resolve(r.result.split(',')[1]);r.onerror=reject;r.readAsDataURL(file)});
  return {population:population(),dataSource:sourceValue(),directory:$('directory').value.trim(),label:$('label').value.trim(),filename:file.name,logDirectory:$('logDirectory').value.trim(),content};
}
$('file').addEventListener('change',()=>{const name=$('file').files[0]?.name;const stem=name?.replace(/\.[^.]+$/,'')||'';$('fileLabel').textContent=name || '选择表型文件';$('label').value=stem.replace(/_Re-fitting$/i,'');$('batchLabel').textContent='本次批次：'+(stem||'待选择表型文件');});
$('directory').addEventListener('input',()=>{$('remotePath').textContent='/data9/home/yzhao/GWAS_IRGSP1.0/'+($('directory').value||'…')});
$('preview').onclick=async()=>{try{const p=await api('/api/preview',await payload());$('previewInfo').textContent=`群体：${p.population} · 标签：${p.label}\n预计 ${p.traits.length} 个性状 / ${p.traits.length*12} 个GWAS作业 / ${p.traits.length} 个绘图作业\n数据源：${JSON.stringify(p.source,null,2)}\n远程目录：${p.remotePath}\n本地日志：${p.logPath}\n性状：\n${p.traits.join('\n')}\n（最终以服务器实际生成并校验的TPED清单为准）`;$('previewText').textContent=p.workflow;$('previewDialog').showModal()}catch(e){message(e.message,true)}};
$('closePreview').onclick=()=>$('previewDialog').close();
$('form').onsubmit=async event=>{event.preventDefault();$('submit').disabled=true;try{const data=await payload();message('正在上传并启动 SecureCRT…');const run=await api('/api/run',data);selected=run.id;message('已启动连接；若需要密码或主机确认，请在 SecureCRT 中完成。');await refresh()}catch(e){message(e.message,true)}finally{$('submit').disabled=false}};
$('history').onchange=()=>{selected=$('history').value;render()};
function render(){const run=runs.find(r=>r.id===selected);if(!run)return;$('statusBadge').textContent=labels[run.status]||run.status;$('runInfo').textContent=`${run.population} · ${run.filename}\n${run.directory}\n${run.detail}\n日志：${run.logPath}`;$('runInfo').style.whiteSpace='pre-wrap';$('log').textContent=run.log||'SecureCRT 正在连接。身份验证如需交互，请打开 SecureCRT 完成。';$('download').disabled=false;
  let stage=run.detail==='prepare'?0:run.detail==='admin'?((run.log.includes('phase=plot')||run.log.includes('PHASE: plot'))?2:1):run.detail==='final'?3:run.status==='complete'?5:-1;
  [...$('steps').children].forEach((li,i)=>{li.className=stage>i?'done':stage===i?'active':''});
}
async function refresh(){try{runs=await api('/api/runs');if(!selected&&runs.length)selected=runs[0].id;$('history').replaceChildren();if(!runs.length){const o=new Option('暂无任务','');$('history').add(o)}for(const r of runs){$('history').add(new Option(`${r.filename} · ${labels[r.status]||r.status}`,r.id));const previous=lastStatus.get(r.id);if(previous&&previous!==r.status&&['complete','attention','failed'].includes(r.status)){message(r.status==='complete'?'分析结束，日志已保存到 '+r.logPath:'任务需要处理，请查看运行日志和 SecureCRT。',r.status!=='complete');document.title=(r.status==='complete'?'✓ 分析完成':'⚠ 需要处理')+' · GWAS'}lastStatus.set(r.id,r.status)}$('history').value=selected;render()}catch(e){$('connection').textContent='本机后台连接中断；请检查服务，勿重复提交。'}}
$('download').onclick=()=>{const r=runs.find(r=>r.id===selected);if(!r)return;const a=document.createElement('a');a.href=URL.createObjectURL(new Blob([r.header+'\n'+r.log+'\n状态：'+r.status+' '+r.detail],{type:'text/plain;charset=utf-8'}));a.download=r.filename.replace(/\.[^.]+$/,'')+'_GWASrun.log';a.click();setTimeout(()=>URL.revokeObjectURL(a.href),1000)};
api('/api/config').then(c=>{config=c;let saved={};try{saved=JSON.parse(localStorage.getItem('gwas-preferences')||'{}')}catch{}$('logDirectory').value=saved.logDirectory||c.defaultLogDirectory;$('directory').value=saved.directory||'';const radio=[...document.querySelectorAll('input[name=population]')].find(r=>r.value===saved.population);if(radio)radio.checked=true;if(c.dataSources)showSource();$('directory').dispatchEvent(new Event('input'));$('bookmarkUrl').href=c.bookmarkUrl;$('bookmarkUrl').textContent=c.bookmarkUrl;$('connection').textContent=c.securecrt?'SecureCRT 已就绪 · 保存会话 admin2':'未找到 SecureCRT，请检查安装路径';$('submit').disabled=!c.securecrt||c.version<4;if(c.version<4){$('preview').disabled=true;message('后台升级中，请稍后刷新页面。')}refresh();setInterval(refresh,5000)}).catch(e=>message('无法连接本机后台：'+e.message,true));
for(const element of [$('directory'),$('logDirectory'),...document.querySelectorAll('input[name=population]')]){element.addEventListener('change',()=>{try{localStorage.setItem('gwas-preferences',JSON.stringify({directory:$('directory').value,logDirectory:$('logDirectory').value,population:document.querySelector('input[name=population]:checked').value}))}catch{}})}
