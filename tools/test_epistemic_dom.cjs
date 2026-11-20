const fs=require('fs');const assert=require('assert');const{JSDOM,VirtualConsole}=require('jsdom');
let checks=[];for(const lang of ['es','en']){
let errors=[];const vc=new VirtualConsole();vc.on('jsdomError',e=>errors.push(e.message));vc.on('error',e=>errors.push(String(e)));
const dom=new JSDOM(fs.readFileSync(require('path').join(__dirname,'../ui/epistemic/epistemic_graph_'+lang+'.html'),'utf8'),{runScripts:'dangerously',pretendToBeVisual:true,url:'https://example.org/',virtualConsole:vc,beforeParse(w){w.matchMedia=()=>({matches:false});w.HTMLElement.prototype.scrollIntoView=()=>{};w.SVGElement.prototype.setPointerCapture=()=>{};w.SVGElement.prototype.hasPointerCapture=()=>false;}});
const w=dom.window,d=w.document;
assert.equal(d.documentElement.dataset.ready,'true',errors.join('\n'));
assert(!d.querySelector('#atlas').hidden);assert(d.body.classList.contains('workflow'));assert.equal(d.querySelectorAll('.node').length,36);d.querySelector('#overviewTab').click();assert(!d.querySelector('#overview').hidden);assert(d.querySelector('#hero-title').textContent.length>30);
d.querySelector('#start').click();assert(!d.querySelector('#story').hidden);assert(d.querySelector('#story-question').textContent.length>20);
for(let i=0;i<6;i++){w.eval(`openChapter(${i}, false)`);assert(d.querySelector('#inspector h2').textContent);assert.equal(d.querySelectorAll('.node').length,6);}
d.querySelector('#pause').click();assert(w.eval('playing'));d.querySelector('#pause').click();assert(!w.eval('playing'));
d.querySelector('#openAtlas').click();assert.equal(d.querySelectorAll('.node').length,36);
d.querySelector('#search').value='Cassini';d.querySelector('#search').dispatchEvent(new w.Event('input'));assert(d.querySelectorAll('.node').length<36);
d.querySelector('#language').value=lang==='es'?'en':'es';d.querySelector('#language').dispatchEvent(new w.Event('change'));assert.equal(d.documentElement.lang,lang==='es'?'en':'es');
d.querySelector('#help').click();assert(!d.querySelector('#modal').hidden);d.querySelector('#closeModal').click();assert(d.querySelector('#modal').hidden);
d.querySelector('#workflowTab').click();w.eval("select('gp_result')");assert.equal(d.querySelectorAll('.node.selected').length,1);assert(d.querySelectorAll('.edge.relevant').length>0);assert(d.querySelectorAll('.node.dim').length>0);
w.eval("dragStart({button:0,pointerType:'mouse',pointerId:1,clientX:200,clientY:200,stopPropagation(){}},'gp_result'); move({pointerId:1,clientX:260,clientY:240});");assert(w.eval('positions.gp_result.x')>0);w.eval("up({pointerId:1})");d.querySelector('#reset').click();assert.equal(Object.keys(w.eval('positions')).length,0);assert(d.querySelector('#minimap').children.length>36);
assert.deepEqual(errors,[]);checks.push({lang,initialization:true,six_chapters:true,play_pause:true,search:true,language:true,modal:true,workflow:true,drag:true,connections_highlighted:true,minimap:true});dom.window.close();
const staticDom=new JSDOM(fs.readFileSync(require('path').join(__dirname,'../ui/epistemic/epistemic_graph_'+lang+'.html'),'utf8'));assert(staticDom.window.document.querySelector('#overview').textContent.length>2000);staticDom.window.close();}
console.log(JSON.stringify({scope:'DOM unit tests; not browser rendering or mobile-device QA',checks},null,2));
