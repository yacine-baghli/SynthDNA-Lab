let selectedProfile='twist',selectedSize=100000,selectedFormat='fasta',selectedSeqBackend='illumina',profiles={},techClasses={},seqBackends={},currentJobId=null;
document.addEventListener('DOMContentLoaded',()=>{loadProfiles();loadSeqBackends();initSliders();initSizeButtons();initFormatButtons();initDNABackground()});

async function loadProfiles(){try{const r=await fetch('/api/profiles');const d=await r.json();profiles=d.profiles;techClasses=d.tech_classes;renderProfileCards();renderProfileDetails();renderCompareSelects();drawErrorCharts()}catch(e){console.error(e)}}
async function loadSeqBackends(){try{const r=await fetch('/api/sequencing_backends');seqBackends=await r.json();renderSeqBackends()}catch(e){console.error(e)}}

function renderProfileCards(){const c=document.getElementById('profile-groups');c.innerHTML='';
const groups={};for(const[k,p]of Object.entries(profiles)){const tc=p.tech_class||'chemical';if(!groups[tc])groups[tc]=[];groups[tc].push({key:k,...p})}
for(const[tc,items]of Object.entries(groups)){const cls=techClasses[tc]||{label:tc,color:'#888'};
const g=document.createElement('div');g.className='tech-group';
g.innerHTML=`<div class="tech-group-label"><span class="tech-dot" style="background:${cls.color}"></span>${cls.label}</div><div class="tech-group-cards" id="tg-${tc}"></div>`;
c.appendChild(g);const row=g.querySelector('.tech-group-cards');
items.forEach(p=>{const card=document.createElement('div');card.className=`profile-card ${p.key===selectedProfile?'active':''}`;card.dataset.key=p.key;
card.style.borderColor=p.key===selectedProfile?cls.color:'';
card.innerHTML=`<div class="card-name">${p.name.split('(')[0].split('/')[0].trim()}</div><div class="card-platform">${(p.platform||'').split(' ').slice(0,2).join(' ')}</div>`;
card.addEventListener('click',()=>selectProfile(p.key,cls.color));row.appendChild(card)})}
updateProfileDesc()}

function selectProfile(k,color){selectedProfile=k;document.querySelectorAll('.profile-card').forEach(c=>{c.classList.remove('active');c.style.borderColor=''});
const el=document.querySelector(`.profile-card[data-key="${k}"]`);if(el){el.classList.add('active');el.style.borderColor=color||'var(--accent-green)'}
updateProfileDesc();renderProfileDetails();drawErrorCharts()}

function updateProfileDesc(){const d=document.getElementById('profile-desc');const p=profiles[selectedProfile];
if(p){d.textContent=p.description;const tc=techClasses[p.tech_class];if(tc)d.style.borderLeftColor=tc.color}}

function renderSeqBackends(){const row=document.getElementById('seq-backend-row');row.innerHTML='';
for(const[k,b]of Object.entries(seqBackends)){const btn=document.createElement('button');btn.className=`seq-btn ${k===selectedSeqBackend?'active':''}`;btn.dataset.key=k;btn.textContent=b.name;
btn.addEventListener('click',()=>{selectedSeqBackend=k;document.querySelectorAll('.seq-btn').forEach(x=>x.classList.remove('active'));btn.classList.add('active');updateSeqDesc()});row.appendChild(btn)}
updateSeqDesc()}

function updateSeqDesc(){const d=document.getElementById('seq-backend-desc');const b=seqBackends[selectedSeqBackend];if(b)d.textContent=b.description}

function renderProfileDetails(){const grid=document.getElementById('profile-detail-grid');const p=profiles[selectedProfile];if(!p)return;
grid.innerHTML=`<div class="detail-card"><h4>Synthesis Errors (Layer 1)</h4>${dr('Deletion rate',sn(p.synth_p_del)+'/base')}${dr('Substitution rate',sn(p.synth_p_sub)+'/base')}${dr('Insertion rate',sn(p.synth_p_ins)+'/base')}${dr('Position slope',p.synth_position_slope.toFixed(2))}${dr('Homopoly alpha',p.homopoly_alpha.toFixed(2))}${dr('Homopoly beta',p.homopoly_beta.toFixed(2))}</div>
<div class="detail-card"><h4>PCR Errors (Layer 2)</h4>${dr('Error/base',sn(p.pcr_error_per_base))}${dr('Sub fraction',(p.pcr_sub_frac*100).toFixed(0)+'%')}${dr('Del fraction',(p.pcr_del_frac*100).toFixed(0)+'%')}${dr('Ins fraction',(p.pcr_ins_frac*100).toFixed(0)+'%')}${dr('Cycles',p.pcr_cycles)}</div>
<div class="detail-card"><h4>Sequencing (Layer 3)</h4>${dr('Sub rate',sn(p.seq_p_sub)+'/base')}${dr('Pos decay',p.seq_position_decay.toFixed(2))}${dr('Trans bias',(p.seq_transition_bias*100).toFixed(0)+'%')}</div>
<div class="detail-card"><h4>Constraints</h4>${dr('GC range',(p.gc_range[0]*100).toFixed(0)+'-'+(p.gc_range[1]*100).toFixed(0)+'%')}${dr('GC optimal',(p.gc_optimal[0]*100).toFixed(0)+'-'+(p.gc_optimal[1]*100).toFixed(0)+'%')}${dr('Max homopoly',p.max_homopolymer+' bp')}${dr('Max dinuc',p.max_dinuc_repeat)}${dr('Oligo range',p.min_oligo_len+'-'+p.max_oligo_len+' nt')}</div>`}

function dr(k,v){return`<div class="detail-row"><span class="detail-key">${k}</span><span class="detail-val">${v}</span></div>`}
function sn(n){if(n===0)return'0';const e=Math.floor(Math.log10(Math.abs(n)));const m=n/Math.pow(10,e);return m.toFixed(1)+'e'+e}

function initSliders(){const ls=document.getElementById('target-len'),lb=document.getElementById('len-badge');ls.addEventListener('input',()=>lb.textContent=ls.value+' nt');
const ws=document.getElementById('num-workers'),wb=document.getElementById('workers-badge');ws.addEventListener('input',()=>wb.textContent=ws.value)}

function initSizeButtons(){document.querySelectorAll('.size-btn').forEach(b=>b.addEventListener('click',()=>{document.querySelectorAll('.size-btn').forEach(x=>x.classList.remove('active'));b.classList.add('active');selectedSize=parseInt(b.dataset.size);document.getElementById('custom-size').value=''}));
document.getElementById('custom-size').addEventListener('input',e=>{const v=parseInt(e.target.value);if(v&&v>0){document.querySelectorAll('.size-btn').forEach(x=>x.classList.remove('active'));selectedSize=v}})}

function initFormatButtons(){document.querySelectorAll('.format-btn').forEach(b=>b.addEventListener('click',()=>{document.querySelectorAll('.format-btn').forEach(x=>x.classList.remove('active'));b.classList.add('active');selectedFormat=b.dataset.format}))}

function switchTab(t){document.querySelectorAll('.tab').forEach(x=>x.classList.remove('active'));document.querySelectorAll('.tab-content').forEach(x=>x.classList.remove('active'));document.querySelector(`.tab[data-tab="${t}"]`).classList.add('active');document.getElementById('tab-'+t).classList.add('active')}

async function previewSample(){const btn=document.getElementById('btn-preview');btn.disabled=true;btn.innerHTML='<span class="spinner"></span>';switchTab('preview');
try{const r=await fetch('/api/preview',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({profile:selectedProfile,sequencing_backend:selectedSeqBackend,target_len:parseInt(document.getElementById('target-len').value)})});renderPreview(await r.json())}catch(e){console.error(e)}finally{btn.disabled=false;btn.innerHTML='<span class="btn-icon">&#128300;</span> Preview'}}

function renderPreview(d){document.getElementById('empty-state').classList.add('hidden');document.getElementById('preview-content').classList.remove('hidden');
document.getElementById('preview-stats').innerHTML=`<div class="stat-card"><div class="stat-value">${d.center_len}</div><div class="stat-label">Length</div></div><div class="stat-card"><div class="stat-value">${(d.gc_content*100).toFixed(1)}%</div><div class="stat-label">GC</div></div><div class="stat-card"><div class="stat-value">${d.max_homopolymer}</div><div class="stat-label">Max HP</div></div><div class="stat-card"><div class="stat-value">${d.num_traces}</div><div class="stat-label">Traces</div></div>`;
document.getElementById('center-meta').textContent=d.center_len+' nt';document.getElementById('center-seq').innerHTML=colorDNA(d.center);
const c=document.getElementById('traces-container');c.innerHTML='';
d.traces.forEach((t,i)=>{const diff=t.len_diff;const ds=diff>0?`-${diff} del`:diff<0?`+${Math.abs(diff)} ins`:'exact';
const b=document.createElement('div');b.className='sequence-block';b.innerHTML=`<div class="seq-header"><span class="seq-label trace-label">TRACE ${i+1}</span><span class="seq-meta">${t.length} nt | ${ds} | GC ${(t.gc_content*100).toFixed(1)}%</span></div><pre class="seq-display">${colorDNA(t.sequence)}</pre>`;c.appendChild(b)})}

function colorDNA(s){return s.split('').map(b=>`<span class="base-${b}">${b}</span>`).join('')}

async function startGeneration(){const btn=document.getElementById('btn-generate');btn.disabled=true;btn.innerHTML='<span class="spinner"></span>';switchTab('generation');
document.getElementById('gen-empty-state').classList.add('hidden');document.getElementById('gen-content').classList.remove('hidden');document.getElementById('gen-results').classList.add('hidden');
document.getElementById('progress-fill').style.width='50%';document.getElementById('progress-text').textContent='...';document.getElementById('gen-message').textContent='Generating samples...';
try{const r=await fetch('/api/generate',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({profile:selectedProfile,sequencing_backend:selectedSeqBackend,target_len:parseInt(document.getElementById('target-len').value),dataset_size:selectedSize,output_format:selectedFormat,min_traces:parseInt(document.getElementById('min-traces').value),max_traces:parseInt(document.getElementById('max-traces').value),num_workers:parseInt(document.getElementById('num-workers').value)})});const d=await r.json();
if(d.status==='complete'||d.fasta_content){renderServerlessResults(d)}
else if(d.job_id){currentJobId=d.job_id;pollJobStatus()}
}catch(e){console.error(e);document.getElementById('gen-message').textContent='Error: '+e.message}finally{btn.disabled=false;btn.innerHTML='<span class="btn-icon">&#9889;</span> Generate'}}

function renderServerlessResults(d){document.getElementById('progress-fill').style.width='100%';document.getElementById('progress-text').textContent='100%';
document.getElementById('gen-message').textContent=d.note||'Complete!';
const c=document.getElementById('gen-results');c.classList.remove('hidden');
let dl='';if(d.fasta_content){const blob=new Blob([d.fasta_content],{type:'text/plain'});const url=URL.createObjectURL(blob);dl=`<a href="${url}" download="synthdna_${selectedProfile}_${d.dataset_size}.fasta" class="download-btn">Download FASTA</a>`}
c.innerHTML=`<div class="result-card"><div class="result-value">${d.dataset_size.toLocaleString()}</div><div class="result-label">Samples</div></div><div class="result-card"><div class="result-value">${d.generation_time_sec}s</div><div class="result-label">Time</div></div><div class="result-card"><div class="result-value">${d.samples_per_sec.toLocaleString()}</div><div class="result-label">Samples/sec</div></div>${dl}`}

async function pollJobStatus(){if(!currentJobId)return;try{const r=await fetch(`/api/status/${currentJobId}`);const d=await r.json();
document.getElementById('progress-fill').style.width=d.progress+'%';document.getElementById('progress-text').textContent=d.progress+'%';document.getElementById('gen-message').textContent=d.message||'...';
if(d.status==='complete'){renderLocalResults(d.result);document.getElementById('btn-generate').disabled=false;document.getElementById('btn-generate').innerHTML='<span class="btn-icon">&#9889;</span> Generate';return}
if(d.status==='error'){document.getElementById('gen-message').textContent='Error: '+(d.error||'Unknown');document.getElementById('btn-generate').disabled=false;document.getElementById('btn-generate').innerHTML='<span class="btn-icon">&#9889;</span> Generate';return}
setTimeout(pollJobStatus,1000)}catch(e){setTimeout(pollJobStatus,2000)}}

function renderLocalResults(r){const c=document.getElementById('gen-results');c.classList.remove('hidden');c.innerHTML=`<div class="result-card"><div class="result-value">${r.dataset_size.toLocaleString()}</div><div class="result-label">Samples</div></div><div class="result-card"><div class="result-value">${r.generation_time_sec}s</div><div class="result-label">Time</div></div><div class="result-card"><div class="result-value">${r.samples_per_sec.toLocaleString()}</div><div class="result-label">Samples/sec</div></div><div class="result-card"><div class="result-value">${r.file_size_mb} MB</div><div class="result-label">File Size</div></div><a href="/api/download/${currentJobId}" class="download-btn">Download ${r.output_format.toUpperCase()}</a>`}

// ══════ COMPARE ══════
function renderCompareSelects(){const a=document.getElementById('compare-a'),b=document.getElementById('compare-b');if(!a||!b)return;a.innerHTML='';b.innerHTML='';
for(const[k,p]of Object.entries(profiles)){a.innerHTML+=`<option value="${k}" ${k==='twist'?'selected':''}>${p.name}</option>`;b.innerHTML+=`<option value="${k}" ${k==='photolitho_ethz'?'selected':''}>${p.name}</option>`}}

async function runComparison(){const a=document.getElementById('compare-a').value,b=document.getElementById('compare-b').value;
const btn=document.querySelector('.btn-compare-go');btn.disabled=true;btn.innerHTML='<span class="spinner"></span>';
try{const r=await fetch('/api/compare',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({profile_a:a,profile_b:b,target_len:parseInt(document.getElementById('target-len').value),n_samples:200})});const d=await r.json();
document.getElementById('compare-results').classList.remove('hidden');drawCompareChart(d,a,b);renderCompareTable(d,a,b)}catch(e){console.error(e)}finally{btn.disabled=false;btn.innerHTML='Compare'}}

function drawCompareChart(data,ka,kb){const canvas=document.getElementById('chart-compare-bars');const ctx=canvas.getContext('2d');const W=canvas.width=canvas.offsetWidth*2;const H=canvas.height=560;ctx.clearRect(0,0,W,H);
const da=data[ka],db=data[kb];const margin={top:50,right:40,bottom:70,left:60};const chartW=W-margin.left-margin.right;const chartH=H-margin.top-margin.bottom;
const metrics=[{label:'Deletion %',a:da.del_frac*100,b:db.del_frac*100},{label:'Substitution %',a:da.sub_frac*100,b:db.sub_frac*100},{label:'Insertion %',a:da.ins_frac*100,b:db.ins_frac*100}];
const maxVal=Math.max(...metrics.map(m=>Math.max(m.a,m.b)))*1.2;const groupW=chartW/metrics.length;const barW=groupW*0.3;
ctx.fillStyle='#8b95a8';ctx.font='600 22px Inter';ctx.textAlign='center';ctx.fillText(`${da.name} vs ${db.name}`,W/2,32);
const tcA=techClasses[da.tech_class]||{color:'#06d6a0'};const tcB=techClasses[db.tech_class]||{color:'#ef476f'};
metrics.forEach((m,i)=>{const gx=margin.left+i*groupW+groupW/2;
const hA=(m.a/maxVal)*chartH;const hB=(m.b/maxVal)*chartH;const xA=gx-barW-4;const xB=gx+4;
const gA=ctx.createLinearGradient(0,margin.top+chartH-hA,0,margin.top+chartH);gA.addColorStop(0,tcA.color);gA.addColorStop(1,tcA.color+'40');ctx.fillStyle=gA;ctx.beginPath();ctx.roundRect(xA,margin.top+chartH-hA,barW,hA,[4,4,0,0]);ctx.fill();
const gB=ctx.createLinearGradient(0,margin.top+chartH-hB,0,margin.top+chartH);gB.addColorStop(0,tcB.color);gB.addColorStop(1,tcB.color+'40');ctx.fillStyle=gB;ctx.beginPath();ctx.roundRect(xB,margin.top+chartH-hB,barW,hB,[4,4,0,0]);ctx.fill();
ctx.fillStyle='#e8ecf4';ctx.font='600 16px "JetBrains Mono"';ctx.fillText(m.a.toFixed(1)+'%',xA+barW/2,margin.top+chartH-hA-8);ctx.fillText(m.b.toFixed(1)+'%',xB+barW/2,margin.top+chartH-hB-8);
ctx.fillStyle='#8b95a8';ctx.font='500 18px Inter';ctx.fillText(m.label,gx,margin.top+chartH+28)});
ctx.fillStyle=tcA.color;ctx.fillRect(margin.left+20,margin.top+chartH+48,14,14);ctx.fillStyle='#e8ecf4';ctx.font='400 15px Inter';ctx.textAlign='left';ctx.fillText(da.name.split('(')[0].trim(),margin.left+42,margin.top+chartH+60);
ctx.fillStyle=tcB.color;ctx.fillRect(margin.left+280,margin.top+chartH+48,14,14);ctx.fillStyle='#e8ecf4';ctx.fillText(db.name.split('(')[0].trim(),margin.left+302,margin.top+chartH+60)}

function renderCompareTable(data,ka,kb){const c=document.getElementById('compare-table');const da=data[ka],db=data[kb];
c.innerHTML=`<div class="compare-card"><h4 style="color:${(techClasses[da.tech_class]||{}).color||'#06d6a0'}">${da.name}</h4>${dr('Total error/nt',sn(da.total_error_rate))}${dr('GC range',(da.gc_range[0]*100).toFixed(0)+'-'+(da.gc_range[1]*100).toFixed(0)+'%')}${dr('Max homopoly',da.max_homopolymer+' bp')}${dr('Max oligo',da.max_oligo_len+' nt')}${dr('Avg len diff',da.avg_len_diff)}${dr('Avg traces',da.avg_traces)}</div>
<div class="compare-card"><h4 style="color:${(techClasses[db.tech_class]||{}).color||'#ef476f'}">${db.name}</h4>${dr('Total error/nt',sn(db.total_error_rate))}${dr('GC range',(db.gc_range[0]*100).toFixed(0)+'-'+(db.gc_range[1]*100).toFixed(0)+'%')}${dr('Max homopoly',db.max_homopolymer+' bp')}${dr('Max oligo',db.max_oligo_len+' nt')}${dr('Avg len diff',db.avg_len_diff)}${dr('Avg traces',db.avg_traces)}</div>`}

// ══════ CHARTS ══════
function drawErrorCharts(){const p=profiles[selectedProfile];if(!p)return;drawErrorDistChart(p);drawPositionBiasChart(p);drawHomopolymerChart(p)}

function drawErrorDistChart(p){const canvas=document.getElementById('chart-error-dist');const ctx=canvas.getContext('2d');const W=canvas.width=canvas.offsetWidth*2;const H=canvas.height=440;ctx.clearRect(0,0,W,H);
const total=p.synth_p_del+p.synth_p_sub+p.synth_p_ins;const bars=[{label:'Deletion',value:p.synth_p_del/total,color:'#ef476f',raw:p.synth_p_del},{label:'Substitution',value:p.synth_p_sub/total,color:'#ffd166',raw:p.synth_p_sub},{label:'Insertion',value:p.synth_p_ins/total,color:'#118ab2',raw:p.synth_p_ins}];
const margin={top:60,right:40,bottom:60,left:60};const chartW=W-margin.left-margin.right;const chartH=H-margin.top-margin.bottom;const barW=chartW/3*0.55;const gap=chartW/3;
const tc=techClasses[p.tech_class]||{color:'#06d6a0',label:''};ctx.fillStyle='#8b95a8';ctx.font='600 22px Inter';ctx.textAlign='center';ctx.fillText('Error Distribution ('+tc.label+')',W/2,32);
bars.forEach((b,i)=>{const x=margin.left+i*gap+(gap-barW)/2;const h=b.value*chartH;const y=margin.top+chartH-h;
const g=ctx.createLinearGradient(x,y,x,y+h);g.addColorStop(0,b.color);g.addColorStop(1,b.color+'60');ctx.fillStyle=g;ctx.beginPath();ctx.roundRect(x,y,barW,h,[6,6,0,0]);ctx.fill();
ctx.fillStyle='#e8ecf4';ctx.font='700 20px "JetBrains Mono"';ctx.textAlign='center';ctx.fillText((b.value*100).toFixed(0)+'%',x+barW/2,y-10);
ctx.fillStyle='#8b95a8';ctx.font='500 18px Inter';ctx.fillText(b.label,x+barW/2,margin.top+chartH+28);
ctx.fillStyle='#5a6478';ctx.font='400 14px "JetBrains Mono"';ctx.fillText(b.raw.toExponential(1)+'/b',x+barW/2,margin.top+chartH+48)})}

function drawPositionBiasChart(p){const canvas=document.getElementById('chart-position-bias');const ctx=canvas.getContext('2d');const W=canvas.width=canvas.offsetWidth*2;const H=canvas.height=440;ctx.clearRect(0,0,W,H);
const margin={top:60,right:40,bottom:50,left:70};const chartW=W-margin.left-margin.right;const chartH=H-margin.top-margin.bottom;
ctx.fillStyle='#8b95a8';ctx.font='600 22px Inter';ctx.textAlign='center';ctx.fillText('Position-Dependent Errors',W/2,32);
const n=50;const synthTotal=p.synth_p_del+p.synth_p_sub+p.synth_p_ins;const maxRate=Math.max(synthTotal*(1+p.synth_position_slope),p.seq_p_sub*(1+p.seq_position_decay))*1.3;
ctx.beginPath();ctx.strokeStyle='#ef476f';ctx.lineWidth=3;for(let i=0;i<n;i++){const f=i/(n-1);const s=1+p.synth_position_slope*(1-f);const r=synthTotal*s;const x=margin.left+f*chartW;const y=margin.top+chartH-(r/maxRate)*chartH;i===0?ctx.moveTo(x,y):ctx.lineTo(x,y)}ctx.stroke();
ctx.beginPath();ctx.strokeStyle='#118ab2';ctx.lineWidth=3;for(let i=0;i<n;i++){const f=i/(n-1);const s=1+p.seq_position_decay*f;const r=p.seq_p_sub*s;const x=margin.left+f*chartW;const y=margin.top+chartH-(r/maxRate)*chartH;i===0?ctx.moveTo(x,y):ctx.lineTo(x,y)}ctx.stroke();
ctx.fillStyle='#5a6478';ctx.font='500 16px Inter';ctx.textAlign='center';ctx.fillText("5' end",margin.left+30,margin.top+chartH+30);ctx.fillText("3' end",margin.left+chartW-30,margin.top+chartH+30);
ctx.fillStyle='#ef476f';ctx.fillRect(margin.left+20,margin.top+10,14,14);ctx.fillStyle='#e8ecf4';ctx.font='400 15px Inter';ctx.textAlign='left';ctx.fillText('Synthesis',margin.left+42,margin.top+22);
ctx.fillStyle='#118ab2';ctx.fillRect(margin.left+140,margin.top+10,14,14);ctx.fillStyle='#e8ecf4';ctx.fillText('Sequencing',margin.left+162,margin.top+22)}

function drawHomopolymerChart(p){const canvas=document.getElementById('chart-homopolymer');const ctx=canvas.getContext('2d');const W=canvas.width=canvas.offsetWidth*2;const H=canvas.height=440;ctx.clearRect(0,0,W,H);
const margin={top:60,right:40,bottom:50,left:70};const chartW=W-margin.left-margin.right;const chartH=H-margin.top-margin.bottom;
ctx.fillStyle='#8b95a8';ctx.font='600 22px Inter';ctx.textAlign='center';ctx.fillText('Homopolymer Scaling',W/2,32);
const a=p.homopoly_alpha,b=p.homopoly_beta,maxK=8;const vals=[];for(let k=1;k<=maxK;k++)vals.push(1+a*(Math.exp(b*Math.max(k-2,0))-1));
const maxVal=Math.max(...vals)*1.15;const barW=chartW/maxK*0.6;const gap=chartW/maxK;
vals.forEach((v,i)=>{const x=margin.left+i*gap+(gap-barW)/2;const h=(v/maxVal)*chartH;const y=margin.top+chartH-h;
const t=Math.min((v-1)/2,1);const r=Math.round(6+t*233),g=Math.round(214-t*143),bl=Math.round(160-t*49);
const gr=ctx.createLinearGradient(x,y,x,y+h);gr.addColorStop(0,`rgb(${r},${g},${bl})`);gr.addColorStop(1,`rgba(${r},${g},${bl},0.3)`);ctx.fillStyle=gr;ctx.beginPath();ctx.roundRect(x,y,barW,h,[6,6,0,0]);ctx.fill();
ctx.fillStyle='#e8ecf4';ctx.font='700 16px "JetBrains Mono"';ctx.textAlign='center';ctx.fillText(v.toFixed(2)+'x',x+barW/2,y-8);
ctx.fillStyle='#8b95a8';ctx.font='500 16px Inter';ctx.fillText((i+1)+'bp',x+barW/2,margin.top+chartH+26)})}

// ══════ DNA BACKGROUND ══════
function initDNABackground(){const canvas=document.getElementById('dna-bg');const ctx=canvas.getContext('2d');let W,H;
function resize(){W=canvas.width=window.innerWidth;H=canvas.height=window.innerHeight}resize();window.addEventListener('resize',resize);
const strands=[];for(let i=0;i<5;i++)strands.push({x:Math.random()*W,speed:0.3+Math.random()*0.4,amp:30+Math.random()*40,freq:0.008+Math.random()*0.006,phase:Math.random()*Math.PI*2,hue:[160,195,340,45,270][i]});
let time=0;function draw(){ctx.clearRect(0,0,W,H);time+=0.01;strands.forEach(s=>{ctx.beginPath();ctx.strokeStyle=`hsla(${s.hue},70%,55%,0.08)`;ctx.lineWidth=1.5;for(let y=0;y<H;y+=3){const x=s.x+Math.sin(y*s.freq+time*s.speed+s.phase)*s.amp;y===0?ctx.moveTo(x,y):ctx.lineTo(x,y)}ctx.stroke();
ctx.beginPath();ctx.strokeStyle=`hsla(${(s.hue+180)%360},70%,55%,0.06)`;for(let y=0;y<H;y+=3){const x=s.x-Math.sin(y*s.freq+time*s.speed+s.phase)*s.amp;y===0?ctx.moveTo(x,y):ctx.lineTo(x,y)}ctx.stroke();
ctx.strokeStyle=`hsla(${s.hue},50%,50%,0.04)`;ctx.lineWidth=1;for(let y=0;y<H;y+=28){const x1=s.x+Math.sin(y*s.freq+time*s.speed+s.phase)*s.amp;const x2=s.x-Math.sin(y*s.freq+time*s.speed+s.phase)*s.amp;ctx.beginPath();ctx.moveTo(x1,y);ctx.lineTo(x2,y);ctx.stroke()}});requestAnimationFrame(draw)}draw()}
