"""Browser behavior checks; no physical validation is inferred."""
from pathlib import Path
from functools import partial
from http.server import ThreadingHTTPServer, SimpleHTTPRequestHandler
from threading import Thread
import json,time
from playwright.sync_api import sync_playwright
R=Path(__file__).resolve().parents[1]
class Quiet(SimpleHTTPRequestHandler):
 def log_message(self,*args):pass
s=ThreadingHTTPServer(('127.0.0.1',0),partial(Quiet,directory=str(R)));Thread(target=s.serve_forever,daemon=True).start()
checks={};errors=[];requests=[]
with sync_playwright() as p:
 b=p.chromium.launch(executable_path='/usr/bin/chromium',args=['--no-sandbox','--disable-dev-shm-usage'])
 page=b.new_page(viewport={'width':1760,'height':1150},device_scale_factor=1)
 page.on('pageerror',lambda e:errors.append(str(e)))
 page.on('request',lambda r:requests.append(r.url))
 page.goto(f'http://127.0.0.1:{s.server_port}/ui/epistemic/epistemic_graph_es.html')
 page.wait_for_function("document.querySelectorAll('.node').length===36")
 checks['initial_36_nodes']=page.locator('.node').count()==36
 page.screenshot(path=str(R/'records/epistemic_overview_es.png'),full_page=True)
 page.locator('[data-id="cassini_check"]').click()
 checks['cassini_failure_visible']='0 / 1024' in page.locator('#inspector').inner_text()
 with page.expect_download() as di:page.locator('[data-artifact="0"]').click()
 download=di.value;download.save_as(str(R/'records/epistemic_test_source.json'))
 checks['embedded_evidence_exact']=(R/'records/epistemic_test_source.json').read_bytes()==(R/'research/revision_v17/results/cassini.json').read_bytes()
 page.select_option('#language','en');checks['english_title']='General Relativity' in page.locator('h1').inner_text()
 page.locator('#tourStart').click();checks['tour_opens']=page.locator('#tour').is_visible()
 page.locator('#next').click();page.locator('#next').click();checks['tour_cassini']='Cassini' in page.locator('#tourTitle').inner_text()
 page.screenshot(path=str(R/'records/epistemic_tour_en.png'),full_page=True)
 page.locator('#play').click();page.wait_for_timeout(650);page.locator('#pause').click();v1=page.evaluate('elapsed');page.wait_for_timeout(650);v2=page.evaluate('elapsed');checks['pause_freezes_timer']=v1==v2 and v1>0
 page.evaluate("elapsed=+$('duration').value-60");page.locator('#pause').click();page.wait_for_timeout(200);checks['automatic_advance']=page.evaluate('tourIndex')==3
 page.locator('#pause').click();page.locator('#closeTour').click();page.locator('[data-route="all"]').click()
 page.locator('#search').fill('Rényi');checks['search_filters']=0<page.locator('.node').count()<36
 page.locator('#search').fill('');page.select_option('#state','unresolved');checks['state_filters']=page.locator('.node').count()==1
 page.select_option('#state','all');page.locator('#fit').click()
 n=page.locator('[data-id="gp_result"]');bb=n.bounding_box();before=n.get_attribute('transform');page.mouse.move(bb['x']+bb['width']/2,bb['y']+bb['height']/2);page.mouse.down();page.mouse.move(bb['x']+bb['width']/2+35,bb['y']+bb['height']/2+25,steps=10);page.mouse.up();checks['card_moves']=page.locator('[data-id="gp_result"]').get_attribute('transform')!=before
 page.locator('#reset').click();checks['layout_resets']=page.locator('[data-id="gp_result"]').get_attribute('transform')==before
 k=page.evaluate('view.k');page.locator('#plus').click();checks['zoom_increases']=page.evaluate('view.k')>k
 page.locator('#depth').click();checks['depth_toggle']=page.locator('#stage').evaluate("e=>e.classList.contains('flat')")
 page.locator('#export').click()
 with page.expect_download() as di:page.locator('#svgExport').click()
 di.value.save_as(str(R/'records/epistemic_export_test.svg'));checks['svg_export']= (R/'records/epistemic_export_test.svg').stat().st_size>1000
 page.locator('#closeModal').click();page.set_viewport_size({'width':390,'height':844});page.locator('#fit').click();page.screenshot(path=str(R/'records/epistemic_mobile_en.png'),full_page=True)
 checks['mobile_no_page_overflow']=page.evaluate('document.documentElement.scrollWidth<=innerWidth')
 page.goto(f'http://127.0.0.1:{s.server_port}/ui/epistemic/epistemic_graph_en.html');checks['english_entrypoint']='General Relativity' in page.locator('h1').inner_text()
 checks['no_external_requests']=all(u.startswith(f'http://127.0.0.1:{s.server_port}') for u in requests)
 checks['no_javascript_errors']=not errors
 b.close()
s.shutdown()
(R/'records/epistemic_browser_tests.json').write_text(json.dumps({'checks':checks,'errors':errors,'scope':'UI behavior only; no native Windows/macOS/iOS execution'},indent=2)+'\n')
print(json.dumps(checks,indent=2));assert all(checks.values()),checks
