"""Optional real-browser checks. Run locally with Playwright + Chromium installed.
This test does not infer physical validation from interface behavior.
"""
from pathlib import Path
import json
from playwright.sync_api import sync_playwright
R=Path(__file__).resolve().parents[1]
checks={};errors=[]
with sync_playwright() as p:
 browser=p.chromium.launch()
 for lang in ['es','en']:
  page=browser.new_page(viewport={'width':1440,'height':1000})
  page.on('pageerror',lambda e:errors.append(str(e)))
  page.goto((R/f'ui/epistemic/epistemic_graph_{lang}.html').as_uri())
  page.wait_for_function("document.documentElement.dataset.ready==='true'")
  checks[lang+'_workflow']=page.locator('#atlas').is_visible()
  page.locator('#overviewTab').click()
  checks[lang+'_overview']=page.locator('#overview').is_visible()
  page.locator('#start').click()
  for i in range(6):
   page.evaluate(f'openChapter({i},false)')
   checks[f'{lang}_chapter_{i}']=page.locator('.node').count()==6 and bool(page.locator('#inspector h2').inner_text())
  page.locator('#pause').click();page.wait_for_timeout(250);page.locator('#pause').click()
  elapsed=page.evaluate('elapsed');page.wait_for_timeout(250)
  checks[lang+'_pause']=elapsed==page.evaluate('elapsed') and elapsed>0
  page.locator('#openAtlas').click();checks[lang+'_36_nodes']=page.locator('.node').count()==36
  page.locator('#search').fill('Rényi');checks[lang+'_filter']=0<page.locator('.node').count()<36
  page.locator('#search').fill('');page.select_option('#route','gp')
  node=page.locator('[data-id="gp_result"]');box=node.bounding_box();before=node.get_attribute('transform')
  page.mouse.move(box['x']+box['width']/2,box['y']+box['height']/2);page.mouse.down();page.mouse.move(box['x']+box['width']/2+30,box['y']+box['height']/2+20,steps=5);page.mouse.up()
  checks[lang+'_drag']=node.get_attribute('transform')!=before
  page.locator('#reset').click();checks[lang+'_reset']=node.get_attribute('transform')==before
  page.set_viewport_size({'width':390,'height':844});page.wait_for_timeout(200)
  checks[lang+'_mobile_width']=page.evaluate('document.documentElement.scrollWidth<=innerWidth')
  page.locator('#mobileMode').click()
  checks[lang+'_mobile_cards']=page.locator('#mobile-path').is_visible()
  page.close()
  offline=browser.new_context(java_script_enabled=False);page=offline.new_page();page.goto((R/f'ui/epistemic/epistemic_graph_{lang}.html').as_uri())
  checks[lang+'_without_js']=len(page.locator('#overview').inner_text())>2000;offline.close()
 browser.close()
result={'scope':'UI behavior only','checks':checks,'errors':errors}
(R/'records/epistemic_browser_checks.json').write_text(json.dumps(result,indent=2)+'\n')
assert all(checks.values()) and not errors,result
print(json.dumps(result,indent=2))
