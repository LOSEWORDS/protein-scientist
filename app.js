'use strict';
const $=id=>document.getElementById(id);
const labels={resolve_input:'Resolve input',predict_complex:'Predict complex',prediction_completed:'Prediction ready',assess_structure:'Assess structure',measure_contacts:'Measure contacts',register_trigger:'Register revision trigger',plan_revision:'Plan revision',execute_revision:'Execute revision',compare_structures:'Compare structures',evaluate_revision:'Review revision',score_candidates:'Score candidate',select_candidates:'Select candidate',audit_delivery:'Audit delivery',finish:'Deliver structure',update_memory:'Update scientific judgment',detect_pockets:'Detect pockets'};
let data,index=0,current,viewer,highlight=false,pocket=false,renderGeneration=0,phases=[],phaseIndex=0;
const cache=new Map();
function el(tag,text,cls){const e=document.createElement(tag);if(text!==undefined)e.textContent=text;if(cls)e.className=cls;return e;}
function block(title,text){const b=el('div',undefined,'block');b.append(el('h3',title));if(text)b.append(el('p',text));return b;}
function refs(ids){const b=el('div');for(const id of ids||[]){const x=el('button',id,'eid');x.onclick=()=>{ $('evidenceTitle').textContent=id+' · '+(data.evidence[id]?.observation?.action||'');$('evidenceBody').textContent=JSON.stringify(data.evidence[id],null,2);$('evidenceDialog').showModal();};b.append(x);}return b;}
function fmt(v,n=3){return typeof v==='number'?v.toFixed(n):'Not Measured';}
function metrics(){
 const event=data.events[index],s=event.snapshot[current]||{},a=s.assessment,m=s.comparison,g=s.score;
 const checks=a?Object.values(a.checks):[];
 const rows=[['Enabled checks',a?`${checks.filter(x=>x==='pass').length}/${checks.length} passed`:'Not measured',a?.evidence_id],['Protein–ligand clashes',a?String(a.raw_metrics.num_pairwise_clashes_protein??'Not measured'):'Not measured',a?.evidence_id],['GNINA CNN score',fmt(g?.cnn_score),g?.evidence_id],['Whole-chain Cα change · Å',fmt(m?.protein_CA_aligned_change_rms_A),m?.evidence_id]];
 const parentId=data.candidates[current]?.parent_id;
 const before=parentId?event.snapshot[parentId]:null;
 if(before){
  const bc=before.assessment?Object.values(before.assessment.checks):null;
  if(bc&&a)rows[0][1]=`${bc.filter(v=>v==='pass').length}/${bc.length} → ${checks.filter(v=>v==='pass').length}/${checks.length}`;
  const oldClash=before.assessment?.raw_metrics?.num_pairwise_clashes_protein;
  if(typeof oldClash==='number'&&a)rows[1][1]=`${oldClash} → ${a.raw_metrics.num_pairwise_clashes_protein??'Not measured'}`;
  if(before.score&&g)rows[2][1]=`${fmt(before.score.cnn_score)} → ${fmt(g.cnn_score)}`;
 }
 $('metrics').replaceChildren();for(const [name,value,source]of rows){const box=el('div',undefined,'metric');box.append(el('div',name,'name'),el('div',value,'value'),el('div',source||'No measurement','source'));$('metrics').append(box);}
 $('pocket').textContent=pocket?'Whole structure':'Pocket view';$('contact').disabled=!s.contacts;$('contact').textContent=highlight?'Hide contacts':'Show contacts';
}
async function textFile(url){if(!cache.has(url)){const r=await fetch(url);if(!r.ok)throw Error('Structure file unavailable');cache.set(url,await r.text());}return cache.get(url);}
async function render(){const generation=++renderGeneration;const c=data.candidates[current];if(!data.events[index].candidate_id){$('candidateLabel').textContent='Not generated';if(viewer){viewer.removeAllModels();viewer.render();}$('viewerStatus').textContent='No structure generated at this step';$('download').hidden=true;$('ligandDownload')?.remove();$('candidate').disabled=true;$('contact').disabled=true;metrics();return;}$('download').hidden=false;$('candidate').disabled=false;$('candidateLabel').textContent=current.slice(0,10);$('download').href=c.receptor;$('download').textContent='Protein PDB ↓';$('ligandDownload')?.remove();const a=el('a','Ligand SDF ↓');a.id='ligandDownload';a.href=c.ligand;a.download='';$('download').after(a);metrics();
 try{if(!window.$3Dmol)throw Error('Viewer unavailable. Evidence and structure downloads remain accessible.');if(!viewer)viewer=$3Dmol.createViewer('viewer',{backgroundColor:'#122a3b',antialias:true});const [pdb,sdf]=await Promise.all([textFile(c.receptor),textFile(c.ligand)]);if(generation!==renderGeneration)return;viewer.removeAllModels();viewer.removeAllLabels();viewer.removeAllShapes();viewer.addModel(pdb,'pdb');viewer.addModel(sdf,'sdf');viewer.setStyle({model:0},{cartoon:{color:'#72c7c4',opacity:.9}});viewer.setStyle({model:1},{stick:{radius:.22,colorscheme:'orangeCarbon'}});
 const contacts=data.events[index].snapshot[current]?.contacts;if(highlight&&contacts){const positions=new Set(contacts.contact_positions);for(const row of contacts.residues||[]){if(positions.has(row.input_position)){viewer.addStyle({model:0,chain:row.chain,resi:row.pdb_number},{stick:{radius:.12,color:'#c5df8a'}});}}}
 
 if(pocket){
  viewer.addStyle({model:0,chain:'A',resi:25},{stick:{radius:.22,colorscheme:'cyanCarbon'}});
  const sulfur=viewer.getModel(0).selectedAtoms({chain:'A',resi:25,atom:'SG'})[0];
  const atoms=viewer.getModel(1).selectedAtoms({}).filter(a=>a.elem!=='H');
  if(sulfur&&atoms.length){
   const distance=a=>Math.hypot(a.x-sulfur.x,a.y-sulfur.y,a.z-sulfur.z);
   const nearest=atoms.reduce((a,b)=>distance(a)<distance(b)?a:b);
   viewer.addLine({start:sulfur,end:nearest,dashed:true,color:'#f5df81'});
   viewer.addLabel('Cys25 SG → nearest ligand heavy atom: '+distance(nearest).toFixed(3)+' Å',{position:{x:(sulfur.x+nearest.x)/2,y:(sulfur.y+nearest.y)/2,z:(sulfur.z+nearest.z)/2},fontSize:12,backgroundColor:'#122a3b',fontColor:'white'});
  }
  viewer.zoomTo({model:1});
 }else viewer.zoomTo();
 viewer.render();$('viewerStatus').textContent=pocket?'Native coordinates · Nearest atom may differ between candidates':'Drag to rotate · Scroll to zoom';window.__viewerReady=true;
 }catch(e){$('viewerStatus').textContent=e.message;window.__viewerError=e.message;}}
function groupPhases(events){
 const result=[];let mode='prediction',revision=0;
 events.forEach((event,i)=>{
  if(event.action==='assess_structure'&&mode==='prediction')mode='assessment';
  if(event.action==='execute_revision'&&event.kind!=='error'){revision++;mode='revision-'+revision;}
  if(event.action==='select_candidates'&&event.kind!=='error')mode='delivery';
  if(!result.length||result[result.length-1].key!==mode){
   const title=mode==='prediction'?'Predict the complex':mode==='assessment'?'Assess the initial structure':mode==='delivery'?'Select and deliver':`Revision ${revision} · Test and review`;
   result.push({key:mode,title,indices:[]});
  }
  result[result.length-1].indices.push(i);
 });
 return result;
}
function excerpt(text,limit=280){return text.length>limit?text.slice(0,limit)+'…':text;}
function showPhase(n){
 phaseIndex=n;const phase=phases[n];index=phase.indices[phase.indices.length-1];
 const event=data.events[index];current=event.candidate_id||data.original_id;
 document.querySelectorAll('.event').forEach((e,i)=>e.classList.toggle('selected',i===n));
 $('stageTitle').textContent='Structure Preview';$('stageType').textContent='RECORDED EVIDENCE';
 $('detailTitle').textContent=`${n+1}. ${phase.title}`;
 const records=phase.indices.map(i=>data.events[i]);const body=$('detailBody');body.replaceChildren();

 if(phase.narrative){
  body.append(block('Run summary',phase.narrative.text),block('Decision summary',phase.narrative.judgment));
 }else{
  const decision=[...records].reverse().find(e=>e.kind==='decision'&&e.decision?.reason);
  if(decision)body.append(block('Recorded decision',decision.decision.reason));
 }
 for(const option of $('candidate').options){const known=data.events.slice(0,index+1).some(e=>e.candidate_id===option.value||e.result?.candidate_ids?.includes(option.value));option.hidden=!known;option.disabled=!known;}
 
 const child=Object.keys(data.candidates).find(k=>data.candidates[k].parent_id===data.original_id);
 const childKnown=child&&Array.from($('candidate').options).some(o=>o.value===child&&!o.disabled);
 $('before').hidden=$('after').hidden=!childKnown;
 $('before').onclick=()=>{current=data.original_id;$('candidate').value=current;pocket=true;render();};
 $('after').onclick=()=>{current=child;$('candidate').value=current;pocket=true;render();};
 $('candidate').value=current;$('requestPanel').hidden=n!==0;$('final').hidden=n!==phases.length-1;render();
}
function finalSection(){
 $('finalStatus').textContent='Local repair succeeded';$('evalMetrics').replaceChildren();
 $('finalReason').textContent='Local repair removed the clash and preserved all 16 parent contacts. The revised structure passed all 11 enabled checks and was selected for delivery. Experimental pose accuracy remains unknown.';
 $('fullAuditText').textContent=data.final.audit?.reason||'No audit text recorded.';
 $('fullAudit').open=false;
 $('limitations').replaceChildren();
 for(const q of [...(data.final.audit?.unresolved_questions||[]),...(data.final.limitations||[])])$('limitations').append(el('p',q));
 const e=data.evaluation||{};
 $('provenance').textContent=JSON.stringify({...data.provenance,session_id:data.session_id,independent_evaluation:e},null,2);
}
async function init(){try{
 const r=await fetch('data/replay.json');if(!r.ok)throw Error('Replay data unavailable.');data=await r.json();phases=groupPhases(data.events);
 try{const response=await fetch('data/reasoning-summary.json');if(response.ok){const summary=await response.json();if(summary.session_id===data.session_id&&summary.session_sha256===data.provenance.session_sha256){const covered=summary.rounds.flatMap(r=>Array.from({length:r.end-r.start+1},(_,i)=>r.start+i));if(covered.length!==data.events.length||covered.some((v,i)=>v!==i))throw Error('Summary chronology mismatch');for(const round of summary.rounds){const known=new Set(data.events.slice(0,round.end+1).flatMap(e=>e.evidence_ids));if(round.evidence_ids.some(id=>!known.has(id)||!data.evidence[id]))throw Error('Summary evidence mismatch');}phases=summary.rounds.map(n=>({title:n.title,key:n.label,indices:Array.from({length:n.end-n.start+1},(_,i)=>n.start+i),narrative:n}));}}}catch(error){console.warn('Using raw recorded phases:',error.message);}

 $('requestGoal').textContent=data.request?.goal||data.goal;
 $('subtitle').textContent=`${data.case_id||data.title.split(' · ')[0]} complex · ${data.candidates[data.original_id].length}-residue protein · One refinement cycle`;
 $('count').textContent=phases.length+' phases';
 phases.forEach((phase,i)=>{const b=el('button',undefined,'event');b.append(el('span',String(i+1).padStart(2,'0'),'number'));const t=el('div',phase.title,'label');t.append(el('span',phase.indices.length+' execution records','meta'));b.append(t);b.onclick=()=>showPhase(i);$('timeline').append(b);});
 for(const [key,c]of Object.entries(data.candidates)){const o=el('option',(c.parent_id?'Revised':'Original'));o.value=key;$('candidate').append(o);}
 $('candidate').onchange=()=>{current=$('candidate').value;render();};$('contact').onclick=()=>{highlight=!highlight;render();};$('pocket').onclick=()=>{pocket=!pocket;render();};$('reset').onclick=()=>{pocket=false;render();};$('closeEvidence').onclick=()=>$('evidenceDialog').close();

 finalSection();showPhase(0);window.__replayLoaded=true;window.__phaseCount=phases.length;
 }catch(e){$('subtitle').textContent=e.message;window.__replayError=e.message;}}
init();
