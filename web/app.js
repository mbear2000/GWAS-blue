const $ = id => document.getElementById(id);
let config, runs = [], selected = '', lastStatus = new Map();
const sourceKeys=['version','sample','tped','kinship','maf'];
let savedSources={};
let phenotypeFiles=[];

try { savedSources=JSON.parse(localStorage.getItem('gwas-data-sources')||'{}'); } catch {}

const population=()=>document.querySelector('input[name=population]:checked').value;
const sourceValue=()=>Object.fromEntries(sourceKeys.map(k=>[k,$('source-'+k).value.trim()]));

function message(text,error=false){
  $('message').textContent=text;
  $('message').className=error?'error':'';
}

function showSource(){
  const s=savedSources[population()]||config.dataSources[population()];
  for(const k of sourceKeys)$('source-'+k).value=s[k]||'';
}
for(const k of sourceKeys){
  $('source-'+k).addEventListener('change',()=>{
    savedSources[population()]=sourceValue();
    try{localStorage.setItem('gwas-data-sources',JSON.stringify(savedSources))}catch{}
  });
}
for(const radio of document.querySelectorAll('input[name=population]')){
  radio.addEventListener('change',()=>{
    if(config?.dataSources)showSource();
    updateMode();
  });
}
$('resetSource').onclick=()=>{
  delete savedSources[population()];
  try{localStorage.setItem('gwas-data-sources',JSON.stringify(savedSources))}catch{}
  showSource();
};

const labels={
  connecting:'连接中',running:'运行中',complete:'分析完成',
  attention:'需要处理',failed:'启动失败',confirm:'等待覆盖确认',
  cancelled:'已停止'
};

function injectBlueUI(){
  const fileInput=$('file');
  fileInput.required=false;
  fileInput.multiple=true;

  const box=document.createElement('div');
  box.id='blueExtension';
  box.innerHTML=`
    <div id="phenotypeList" class="blue-file-list"></div>
    <button type="button" class="secondary blue-add" id="addPhenotype">＋ Add phenotype</button>

    <div id="analysisMode" class="blue-mode single">
      分析模式：Single-phenotype GWAS
    </div>

    <div id="bluePanel" class="blue-panel" hidden>
      <div class="blue-title">Multi-environment BLUE + GWAS</div>
      <p>
        系统根据地点和年份自动选择分析模式，不需要手工选择 BLUE 类型。
      </p>
      <div id="blueSummary" class="blue-preview"></div>
    </div>`;

  $('drop').insertAdjacentElement('afterend',box);

  const style=document.createElement('style');
  style.textContent=`
    .blue-add{margin-top:10px;width:100%}
    .blue-mode{margin-top:12px;padding:10px 12px;border-radius:7px;font-size:13px;font-weight:700}
    .blue-mode.single{background:#f2f7f9;color:#49636e}
    .blue-mode.multi{background:#eaf6f3;color:#08796e;border:1px solid #bdddd8}
    .blue-panel{margin-top:12px;padding:15px;border:1px solid #bdddd8;border-radius:8px;background:#f6fbfa}
    .blue-title{font-weight:800;font-size:15px}
    .blue-panel>p{font-size:12px;color:#60757f;line-height:1.6}
    .blue-file-list{display:flex;flex-direction:column;gap:8px;margin-top:12px}
    .blue-file-row{
      display:grid;
      grid-template-columns:minmax(180px,1.35fr) minmax(160px,.8fr) auto;
      gap:8px;align-items:center;padding:9px 10px;
      border:1px solid #dde5e8;border-radius:7px;background:#fff
    }
    .blue-file-row .name{font-size:13px;overflow-wrap:anywhere}
    .blue-file-row input{border:1px solid #cbd8dd;border-radius:5px;padding:8px;font:13px inherit;width:100%;box-sizing:border-box}
    .blue-file-row button{padding:8px 11px;background:#fff;border:1px solid #d8e2e7;color:#8b3a35}
    .blue-preview{white-space:pre-wrap;font:12px/1.65 Consolas,monospace;background:#fff;border-radius:6px;padding:10px;border:1px solid #e0e8eb;margin-top:9px}
    @media(max-width:850px){.blue-file-row{grid-template-columns:1fr}.blue-file-row button{width:100%}}
  `;
  document.head.appendChild(style);
  $('addPhenotype').onclick=()=>fileInput.click();
}

function readAsBase64(file){
  return new Promise((resolve,reject)=>{
    const r=new FileReader();
    r.onload=()=>resolve(r.result.split(',')[1]);
    r.onerror=reject;
    r.readAsDataURL(file);
  });
}
function readFirstLine(file){
  return new Promise((resolve,reject)=>{
    const r=new FileReader();
    r.onload=()=>resolve(String(r.result||'').replace(/^\uFEFF/,'').split(/\r?\n/)[0]||'');
    r.onerror=reject;
    r.readAsText(file.slice(0,65536),'UTF-8');
  });
}
function splitHeader(header){
  const delim=header.includes('\t')?'\t':',';
  return header.split(delim).map(x=>x.replace(/^"|"$/g,'').trim());
}
function inferPrefix(header,filename){
  const trait=(splitHeader(header)[1]||'');
  let m=trait.match(/^([A-Za-z]+20\d{2}(?:_\d+)?_)/);
  if(m)return m[1];
  m=trait.match(/^([A-Za-z]+20\d{2}_)/);
  if(m)return m[1];

  m=filename.match(/(?:^|_)([A-Za-z]+20\d{2}(?:_\d+)?_)/);
  return m?m[1]:'';
}
function cleanedTraits(item){
  return splitHeader(item.header).slice(1).map(x=>{
    if(item.prefix && x.startsWith(item.prefix))return x.slice(item.prefix.length);
    return x;
  });
}
function inferEnvironment(item){
  const trait=(splitHeader(item.header)[1]||'');
  const candidates=[
    (item.prefix||'').replace(/_$/,''),
    trait,
    item.file.name.replace(/\.[^.]+$/,'')
  ];
  const pats=[
    /^(?:GPallInd_|GPallJap_|GPall_)?([A-Za-z]+)(20\d{2})(?:_|$)/,
    /(?:^|_)([A-Za-z]+)(20\d{2})(?:_|$)/
  ];
  for(const s of candidates){
    for(const p of pats){
      const m=s.match(p);
      if(m)return {location:m[1],year:m[2]};
    }
  }
  return null;
}

function automaticPlan(){
  if(phenotypeFiles.length<2)return null;

  const envs=phenotypeFiles.map(inferEnvironment);
  if(envs.some(x=>!x))return null;

  const locs=[...new Map(envs.map(e=>[e.location.toLowerCase(),e.location])).values()]
    .sort((a,b)=>a.localeCompare(b));
  const years=[...new Set(envs.map(e=>e.year))].sort();

  const traits=cleanedTraits(phenotypeFiles[0]);
  const stem=traits.length===1?traits[0]:'';
  const part=stem?stem+'_':'';
  const pop=population();

  const userLabel=$('label').value.trim();
  const fileLabel=userLabel.replaceAll('-','_');

  if(locs.length===1 && years.length>=2){
    return {
      kind:'multi_env',
      title:`${locs[0]} 多年单点 BLUE`,
      label:userLabel,
      output:userLabel?`${pop}_${part}${fileLabel}`:'',
      program:'generate_blue_multi_env_allTrait.R',
      model:'trait ~ Accession + (1 | env)'
    };
  }

  if(locs.length>=2 && years.length===1){
    return {
      kind:'multi_env',
      title:`${years[0]} 单年多点 BLUE`,
      label:userLabel,
      output:userLabel?`${pop}_${part}${fileLabel}`:'',
      program:'generate_blue_multi_env_allTrait.R',
      model:'trait ~ Accession + (1 | env)'
    };
  }

  if(locs.length>=2 && years.length>=2){
    return {
      kind:'locyear',
      title:'多年多点 Loc × Year BLUE',
      label:userLabel,
      output:userLabel?`${pop}_${part}${fileLabel}`:'',
      program:'generate_blue_multiEnvYear_allTrait.R',
      model:'trait ~ Accession + (1 | Loc) + (1 | Year) + (1 | Loc:Year)'
    };
  }

  return null;
}

function renderPhenotypes(){
  const list=$('phenotypeList');
  list.replaceChildren();

  phenotypeFiles.forEach((item,i)=>{
    const row=document.createElement('div');
    row.className='blue-file-row';

    const name=document.createElement('div');
    name.className='name';
    const env=inferEnvironment(item);
    name.textContent=`${i+1}. ${item.file.name}`+(env?`  [${env.location}${env.year}]`:'');

    const prefix=document.createElement('input');
    prefix.value=item.prefix||'';
    prefix.placeholder='环境前缀，如 SJ2023_';
    prefix.title='只删除性状名称开头完全匹配的环境前缀';
    prefix.oninput=()=>{item.prefix=prefix.value;updateMode()};

    const del=document.createElement('button');
    del.type='button';
    del.textContent='删除';
    del.onclick=()=>{phenotypeFiles.splice(i,1);updateMode()};

    row.append(name,prefix,del);
    list.appendChild(row);
  });

  $('fileLabel').textContent=phenotypeFiles.length
    ?`已选择 ${phenotypeFiles.length} 个表型文件`
    :'选择表型文件';
}

function updateMode(){
  renderPhenotypes();

  const multi=phenotypeFiles.length>=2;
  $('bluePanel').hidden=!multi;
  $('analysisMode').className='blue-mode '+(multi?'multi':'single');
  $('analysisMode').textContent='分析模式：'+
    (multi?'Multi-environment BLUE + GWAS':'Single-phenotype GWAS');

  if(phenotypeFiles.length===1){
    const stem=phenotypeFiles[0].file.name.replace(/\.[^.]+$/,'');
    $('label').value=stem.replace(/_Re-fitting$/i,'');
    $('batchLabel').textContent='本次批次：'+stem;
    return;
  }

  if(!multi){
    $('label').value='';
    $('batchLabel').textContent='本次批次：待选择表型文件';
    return;
  }

  const envs=phenotypeFiles.map(inferEnvironment);
  const plan=automaticPlan();
  const locs=[...new Set(envs.filter(Boolean).map(x=>x.location))];
  const years=[...new Set(envs.filter(Boolean).map(x=>x.year))].sort();

  $('batchLabel').textContent=`本次批次：BLUE（${phenotypeFiles.length} 个环境）`;

  const lines=[];
  lines.push(`地点：${locs.join(', ')||'未识别'}`);
  lines.push(`年份：${years.join(', ')||'未识别'}`);
  lines.push(`清理后性状数：${phenotypeFiles.length?cleanedTraits(phenotypeFiles[0]).length:0}`);
  lines.push('');

  phenotypeFiles.forEach((x,i)=>{
    const e=envs[i];
    lines.push(`${x.file.name}`);
    lines.push(`  环境：${e?e.location+e.year:'无法识别'}`);
    lines.push(`  去除前缀：${x.prefix||'(无)'}`);
    lines.push(`  性状示例：${cleanedTraits(x).slice(0,3).join(', ')}`);
  });

  lines.push('');
  if(plan){
    lines.push(`自动分析：${plan.title}`);
    lines.push(`使用程序：${plan.program}`);
    lines.push(`本次标签：${plan.label||'(请在上方填写)'}`);
    lines.push(`BLUE输出：${plan.output?plan.output+'.txt / .csv':'填写本次标签后生成'}`);
    lines.push(`模型：${plan.model}`);
  }else{
    lines.push('当前环境组合无法确定 BLUE 分析模式。');
  }

  $('blueSummary').textContent=lines.join('\n');
}

injectBlueUI();

function injectInflationStep(){
  const steps=$('steps');
  if(!steps || document.getElementById('inflationStep'))return;

  const complete=[...steps.children].find(
    li=>li.dataset.stage==='complete'
  );
  if(!complete)return;

  const li=document.createElement('li');
  li.id='inflationStep';
  li.dataset.stage='inflation';
  li.innerHTML=`
    <b>5</b>
    <div>
      <strong>计算基因组膨胀系数</strong>
      <span>读取归档 output/*.ps · 计算 Lambda GC 并保存 inflationFactor</span>
    </div>
  `;
  steps.insertBefore(li,complete);
}
injectInflationStep();

function injectRunActionButtons(){
  if(document.getElementById('inflationOnly'))return;

  const box=document.createElement('div');
  box.id='runActionButtons';
  box.style.display='flex';
  box.style.flexWrap='wrap';
  box.style.gap='8px';
  box.style.marginTop='12px';
  box.innerHTML=`
    <button type="button" id="inflationOnly" class="secondary" hidden>
      只计算膨胀系数（使用已有 .ps）
    </button>
    <button type="button" id="cancelLocal" class="secondary" hidden>
      停止本地监控 / 清除运行状态
    </button>
  `;

  const recovery=$('recoveryHelp');
  const host=recovery?.parentElement || $('runInfo')?.parentElement || document.body;
  host.insertBefore(box,recovery||null);

  $('inflationOnly').onclick=async()=>{
    const r=runs.find(x=>x.id===selected);
    if(!r)return;
    if(!confirm(
      '将直接读取该任务已归档的 EMMAX output/*.ps 计算膨胀系数，'+
      '不会重新运行 prepare、qsub、GWAS、绘图或 peak SNP。继续吗？'
    ))return;

    $('inflationOnly').disabled=true;
    try{
      const created=await api('/api/inflation-only',{id:r.id});
      selected=created.id;
      message('已启动“只计算膨胀系数”；正在连接 fat2 并复用已有 .ps。');
      await refresh();
    }catch(e){
      message(e.message,true);
    }finally{
      $('inflationOnly').disabled=false;
    }
  };

  $('cancelLocal').onclick=async()=>{
    const r=runs.find(x=>x.id===selected);
    if(!r)return;
    if(!confirm(
      '这只会停止/清除 Windows 网页中的本地任务状态。'+
      '不会自动 qdel 服务器队列，也不会删除服务器结果。'+
      '请先确认远程任务确实已经中断或不再需要监控。继续吗？'
    ))return;

    $('cancelLocal').disabled=true;
    try{
      await api('/api/cancel-local',{id:r.id});
      message('本地网页任务状态已标记为“已停止”。');
      await refresh();
    }catch(e){
      message(e.message,true);
    }finally{
      $('cancelLocal').disabled=false;
    }
  };
}
injectRunActionButtons();



function injectOverwriteUI(){
  if(document.getElementById('overwriteBox'))return;
  const box=document.createElement('div');
  box.id='overwriteBox';
  box.hidden=true;
  box.style.marginTop='12px';
  box.style.padding='14px';
  box.style.border='1px solid #e6b84f';
  box.style.borderRadius='8px';
  box.style.background='#fff8e6';
  box.innerHTML=`
    <div id="overwriteText" style="margin-bottom:10px;line-height:1.6"></div>
    <button type="button" id="overwriteProceed">继续并覆盖</button>
    <button type="button" id="overwriteStop" style="margin-left:10px">停止本次分析</button>
  `;
  const recovery=document.getElementById('recoveryHelp');
  const host=recovery?.parentElement || document.getElementById('runInfo')?.parentElement || document.body;
  host.insertBefore(box,recovery||null);

  document.getElementById('overwriteProceed').onclick=async()=>{
    const r=runs.find(x=>x.id===selected);
    if(!r)return;
    const btn=document.getElementById('overwriteProceed');
    btn.disabled=true;
    try{
      await api('/api/overwrite-decision',{id:r.id,proceed:true});
      message('已确认覆盖；正在重新连接并继续。只覆盖本次同名目标。');
      await refresh();
    }catch(e){
      message(e.message,true);
    }finally{
      btn.disabled=false;
    }
  };

  document.getElementById('overwriteStop').onclick=async()=>{
    const r=runs.find(x=>x.id===selected);
    if(!r)return;
    const btn=document.getElementById('overwriteStop');
    btn.disabled=true;
    try{
      await api('/api/overwrite-decision',{id:r.id,proceed:false});
      message('已停止本次分析；未进入 prepare/qsub。');
      await refresh();
    }catch(e){
      message(e.message,true);
    }finally{
      btn.disabled=false;
    }
  };
}
injectOverwriteUI();


$('label').addEventListener('input',()=>{ if(phenotypeFiles.length>=2) updateMode(); });

$('file').addEventListener('change',async()=>{
  const files=[...$('file').files];
  $('file').value='';

  for(const file of files){
    if(file.size>20*1024*1024){
      message(`文件不能超过20 MB：${file.name}`,true);
      continue;
    }
    if(phenotypeFiles.some(x=>x.file.name===file.name)){
      message(`已添加同名文件：${file.name}`,true);
      continue;
    }
    try{
      const header=await readFirstLine(file);
      phenotypeFiles.push({
        file,
        header,
        prefix:inferPrefix(header,file.name)
      });
    }catch(e){
      message(`读取文件失败：${file.name} · ${e.message}`,true);
    }
  }
  updateMode();
});

async function payload(){
  if(!$('form').reportValidity())throw Error('请补全分析配置');
  if(!phenotypeFiles.length)throw Error('请至少选择一个表型文件');

  const base={
    population:population(),
    dataSource:sourceValue(),
    directory:$('directory').value.trim(),
    logDirectory:$('logDirectory').value.trim()
  };

  if(phenotypeFiles.length===1){
    const item=phenotypeFiles[0];
    const content=await readAsBase64(item.file);
    return {
      ...base,
      label:$('label').value.trim(),
      filename:item.file.name,
      content,
      phenotypes:[{filename:item.file.name,content,prefix:''}]
    };
  }

  const userLabel=$('label').value.trim();
  if(!userLabel)throw Error('请填写“本次标签”');

  const plan=automaticPlan();
  if(!plan)throw Error('当前文件组合无法确定 BLUE 分析模式，请检查地点/年份识别。');

  const phenotypes=[];
  for(const item of phenotypeFiles){
    phenotypes.push({
      filename:item.file.name,
      prefix:item.prefix||'',
      content:await readAsBase64(item.file)
    });
  }

  return {
    ...base,
    label:userLabel,
    phenotypes
  };
}

$('directory').addEventListener('input',()=>{
  $('remotePath').textContent='/data9/home/yzhao/GWAS_IRGSP1.0/'+($('directory').value||'…');
});

async function api(url,data){
  const response=await fetch(url,data?{
    method:'POST',
    headers:{'Content-Type':'application/json','X-GWAS-Token':config.token},
    body:JSON.stringify(data)
  }:{});
  const result=await response.json();
  if(!response.ok)throw Error(result.error||'请求失败');
  return result;
}

$('preview').onclick=async()=>{
  try{
    const data=await payload();
    const p=await api('/api/preview',data);

    let extra='';
    if(phenotypeFiles.length>=2){
      const plan=automaticPlan();
      extra+=`模式：Multi-environment BLUE + GWAS\n`;
      extra+=`自动分析：${plan.title}\n`;
      extra+=`BLUE标签：${plan.label}\n`;
      extra+=`BLUE输出：${plan.output}\n`;
      extra+=`模型：${plan.model}\n`;
    }else{
      extra+='模式：Single-phenotype GWAS\n';
    }

    $('previewInfo').textContent=
      `${extra}群体：${p.population} · 标签：${p.label}\n`+
      `预计 ${p.traits.length} 个GWAS性状 / ${p.traits.length*12} 个GWAS作业 / ${p.traits.length} 个绘图作业\n`+
      `数据源：${JSON.stringify(p.source,null,2)}\n`+
      `远程目录：${p.remotePath}\n`+
      `本地日志：${p.logPath}\n`+
      `最终GWAS性状：\n${p.traits.join('\n')}\n`;

    $('previewText').textContent=p.workflow;
    $('previewDialog').showModal();
  }catch(e){
    message(e.message,true);
  }
};

$('closePreview').onclick=()=>$('previewDialog').close();

$('form').onsubmit=async event=>{
  event.preventDefault();
  $('submit').disabled=true;

  try{
    const data=await payload();
    message(
      phenotypeFiles.length>=2
        ?'正在上传多个表型并启动 BLUE + GWAS…'
        :'正在上传并启动 SecureCRT…'
    );
    const run=await api('/api/run',data);
    selected=run.id;
    message('已启动连接；如需密码或主机确认，请在 SecureCRT 中完成。');
    await refresh();
  }catch(e){
    message(e.message,true);
  }finally{
    $('submit').disabled=false;
  }
};

$('history').onchange=()=>{
  selected=$('history').value;
  render();
};

function simpleRunName(run){
  return run.displayName || (run.filename||'run').replace(/\.[^.]+$/,'');
}

function maybeNotify(run){
  if(!('Notification' in window) || Notification.permission!=='granted')return;
  const name=simpleRunName(run);
  const title=run.status==='complete'?'GWAS 分析完成':'GWAS 任务需要处理';
  const body=run.status==='complete'
    ?`${name} 分析已完成`
    :`${name}：${labels[run.status]||run.status}，请查看日志`;
  try{ new Notification(title,{body}); }catch{}
}

function render(){
  const run=runs.find(r=>r.id===selected);
  if(!run)return;

  $('statusBadge').textContent=labels[run.status]||run.status;
  $('runInfo').textContent=
    `${simpleRunName(run)}\n${run.directory}\n${run.detail}\n日志：${run.logPath}`;
  $('runInfo').style.whiteSpace='pre-wrap';
  $('log').textContent=
    run.log||'SecureCRT 正在连接。身份验证如需交互，请打开 SecureCRT 完成。';
  $('download').disabled=false;
  $('recoveryHelp').hidden=!['attention','failed'].includes(run.status);
  $('recover').hidden=!!run.remoteDirectoryDeleted;

  const inflationOnlyBtn=$('inflationOnly');
  if(inflationOnlyBtn){
    inflationOnlyBtn.hidden =
      !!run.remoteDirectoryDeleted ||
      !['attention','failed','cancelled','complete'].includes(run.status);
  }

  const cancelLocalBtn=$('cancelLocal');
  if(cancelLocalBtn){
    cancelLocalBtn.hidden =
      !['connecting','running','attention','failed'].includes(run.status);
  }

  const overwriteBox=document.getElementById('overwriteBox');
  if(overwriteBox){
    overwriteBox.hidden=run.status!=='confirm';
    const overwriteText=document.getElementById('overwriteText');
    if(overwriteText && run.status==='confirm'){
      overwriteText.textContent=run.detail ||
        '检测到同名输入表型或本次 BLUE/GWAS 输出。确认继续后，只覆盖本次同名目标。';
    }
  }

  let stage=
    run.detail==='prepare'?0:
    run.detail==='admin'?((run.log.includes('phase=plot')||run.log.includes('PHASE: plot'))?2:1):
    run.detail==='inflation'?4:
    run.detail==='final'?(
      run.log.includes('INFLATION_FACTOR_START') ||
      run.log.includes('===== GENOMIC INFLATION FACTOR =====')
        ?4:3
    ):
    run.status==='complete'?6:-1;

  [...$('steps').children].forEach((li,i)=>{
    li.className=stage>i?'done':stage===i?'active':'';
  });
}

$('recover').onclick=async()=>{
  const r=runs.find(r=>r.id===selected);
  if(!r)return;

  $('recover').disabled=true;
  try{
    await api('/api/recover',{id:r.id});
    message('正在重新连接并核实远程阶段；不会重复提交已启动的GWAS。');
    await refresh();
  }catch(e){
    message(e.message,true);
  }finally{
    $('recover').disabled=false;
  }
};

async function refresh(){
  try{
    runs=await api('/api/runs');
    if(!selected&&runs.length)selected=runs[0].id;

    $('history').replaceChildren();
    if(!runs.length)$('history').add(new Option('暂无任务',''));

    for(const r of runs){
      $('history').add(
        new Option(`${simpleRunName(r)} · ${labels[r.status]||r.status}`,r.id)
      );

      const previous=lastStatus.get(r.id);
      if(previous&&previous!==r.status&&['complete','attention','failed'].includes(r.status)){
        message(
          r.status==='complete'
            ?'分析结束，日志已保存到 '+r.logPath
            :'任务需要处理，请查看运行日志和 SecureCRT。',
          r.status!=='complete'
        );
        maybeNotify(r);
      }
      lastStatus.set(r.id,r.status);
    }

    $('history').value=selected;
    render();
  }catch(e){
    $('connection').textContent='本机后台连接中断；请检查服务，勿重复提交。';
  }
}

$('download').onclick=()=>{
  const r=runs.find(r=>r.id===selected);
  if(!r)return;

  const a=document.createElement('a');
  a.href=URL.createObjectURL(
    new Blob(
      [r.header+'\n'+r.log+'\n状态：'+r.status+' '+r.detail],
      {type:'text/plain;charset=utf-8'}
    )
  );

  const base=simpleRunName(r);
  const m=(r.logPath||'').match(/(\d{8}-\d{6})(?:\.log)?$/);
  const timestamp=m?m[1]:'';
  a.download=timestamp ? `${base}_${timestamp}.log` : `${base}.log`;

  document.body.appendChild(a);
  a.click();
  a.remove();
  setTimeout(()=>URL.revokeObjectURL(a.href),1000);
};

api('/api/config').then(c=>{
  config=c;

  let saved={};
  try{
    saved=JSON.parse(localStorage.getItem('gwas-preferences')||'{}');
  }catch{}

  $('logDirectory').value=saved.logDirectory||c.defaultLogDirectory;
  $('directory').value=saved.directory||'';

  const radio=[...document.querySelectorAll('input[name=population]')]
    .find(r=>r.value===saved.population);
  if(radio)radio.checked=true;

  if(c.dataSources)showSource();

  $('directory').dispatchEvent(new Event('input'));
  $('bookmarkUrl').href=c.bookmarkUrl;
  $('bookmarkUrl').textContent=c.bookmarkUrl;

  $('connection').textContent=c.securecrt
    ?'SecureCRT 已就绪 · 保存会话 admin2'
    :'未找到 SecureCRT，请检查安装路径';

  $('submit').disabled=!c.securecrt||c.version<4;
  if(c.version<4){
    $('preview').disabled=true;
    message('后台升级中，请稍后刷新页面。');
  }

  refresh();
  setInterval(refresh,5000);
}).catch(e=>message('无法连接本机后台：'+e.message,true));

for(const element of [
  $('directory'),
  $('logDirectory'),
  ...document.querySelectorAll('input[name=population]')
]){
  element.addEventListener('change',()=>{
    try{
      localStorage.setItem(
        'gwas-preferences',
        JSON.stringify({
          directory:$('directory').value,
          logDirectory:$('logDirectory').value,
          population:document.querySelector('input[name=population]:checked').value
        })
      );
    }catch{}
  });
}
