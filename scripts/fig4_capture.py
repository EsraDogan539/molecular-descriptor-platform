import json
from pathlib import Path
from playwright.sync_api import sync_playwright
out = Path("docs/fig4"); out.mkdir(parents=True, exist_ok=True)
boxes = {}
def box(page, key, text, exact=False):
    loc = page.get_by_text(text, exact=exact)
    try:
        b = loc.first.bounding_box(); boxes[key] = b
    except Exception as e:
        boxes[key] = repr(e)
with sync_playwright() as p:
    br = p.chromium.launch()
    page = br.new_page(viewport={"width": 1440, "height": 2600}, device_scale_factor=2)
    def go(url):
        page.goto(url, timeout=180000)
        w = page.get_by_text("Yes, get this app back up")
        if w.count(): w.first.click()
        page.wait_for_timeout(50000)
    go("https://chalmoldb.streamlit.app/~/+/?page=database")
    for k, t in [("a_caption", "5,029 records ·"), ("a_matching", "matching records"), ("a_replay", "Replay saved query")]:
        box(page, "db_" + k, t)
    page.screenshot(path=str(out / "database.png"), full_page=True)
    (out / "database.txt").write_text(page.locator("body").inner_text())
    go("https://chalmoldb.streamlit.app/~/+/?page=statistics")
    for k, t in [("snapshot", "Release snapshot"), ("distribution", "Distribution overview"), ("composition", "Collection composition"), ("coverage", "Chalcogen coverage")]:
        box(page, "st_" + k, t)
    page.screenshot(path=str(out / "statistics.png"), full_page=True)
    (out / "statistics.txt").write_text(page.locator("body").inner_text())
    go("https://chalmoldb.streamlit.app/~/+/?page=database&record=EXT_0162")
    page.screenshot(path=str(out / "record.png"), full_page=True)
    br.close()
(out / "boxes.json").write_text(json.dumps(boxes, indent=1))
print(json.dumps(boxes, indent=1))
