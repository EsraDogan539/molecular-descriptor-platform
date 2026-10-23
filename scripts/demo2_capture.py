"""Capture v1.2 site states for the demo video (2x DPR, tall viewport) plus element boxes.
Each step uses its own page, closed afterwards, so app sessions do not pile up."""
import json, sys, traceback
from pathlib import Path
from playwright.sync_api import sync_playwright
OUT = Path("docs/demo2"); OUT.mkdir(parents=True, exist_ok=True)
BASE = "https://chalmoldb.streamlit.app/~/+/"
boxes = json.loads((OUT / "boxes.json").read_text()) if (OUT / "boxes.json").exists() else {}
log = []
ONLY = sys.argv[1].split(",") if len(sys.argv) > 1 else None
HIDE = ("header[data-testid='stHeader'],[data-testid='stToolbar'],[data-testid='stDecoration'],#MainMenu,footer,"
        ".stDeployButton,[data-testid='manage-app-button'],[class*='viewerBadge']{display:none!important}")

def bb(page, key, loc):
    try:
        boxes[key] = loc.first.bounding_box(timeout=8000)
    except Exception as e:
        boxes[key] = None; log.append(f"box {key}: {str(e)[:160]}")

def shot(page, name):
    page.screenshot(path=str(OUT / f"{name}.png"))
    (OUT / f"{name}.txt").write_text(page.locator("body").inner_text())
    log.append(f"shot {name}")

def go(page, q, wait=40000):
    for attempt in range(4):
        page.goto(BASE + q, timeout=180000)
        w = page.get_by_text("Yes, get this app back up")
        if w.count(): w.first.click()
        page.wait_for_timeout(wait)
        body = page.locator("body").inner_text()
        if "supposed to happen" not in body:
            break
        log.append(f"ERROR PAGE at {q} (attempt {attempt + 1})")
        page.wait_for_timeout(60000)
    page.add_style_tag(content=HIDE); page.wait_for_timeout(600)

with sync_playwright() as p:
    br = p.chromium.launch()
    ctx = br.new_context(viewport={"width": 1440, "height": 3000}, device_scale_factor=2)

    def step(name, fn):
        if ONLY and name not in ONLY: return
        page = ctx.new_page()
        try:
            fn(page); log.append(f"ok {name}")
        except Exception:
            log.append(f"FAIL {name}: " + traceback.format_exc()[-600:])
            try: shot(page, name + "_fail")
            except Exception: pass
        finally:
            page.close(); ctx.pages and None
            import time; time.sleep(20)

    def about(page):
        go(page, "?page=about"); shot(page, "about")
        bb(page, "about_cite", page.get_by_text("How to cite", exact=False))
    step("about", about)

    def docs(page):
        go(page, "?page=documentation"); shot(page, "docs")
    step("docs", docs)

    def home(page):
        go(page, "?page=home"); shot(page, "home")
        bb(page, "home_explore", page.get_by_text("Explore database", exact=True))
        bb(page, "home_glance", page.get_by_text("Database at a glance"))
        bb(page, "home_polymer_card", page.get_by_text("Polymer collection · 1,323 records"))
        bb(page, "home_cite", page.get_by_text("How to cite ChalMolDB"))
    step("home", home)

    def database(page):
        go(page, "?page=database"); shot(page, "db_0")
        search = page.get_by_placeholder("ID, molecule/system, structure, method, source or reference")
        bb(page, "db_search", search)
        bb(page, "db_count", page.get_by_text("matching records"))
        bb(page, "db_table", page.get_by_test_id("stDataFrame"))
        info = page.evaluate("""() => [...document.querySelectorAll('[data-testid=stSelectbox]')].map(e => {const r=e.getBoundingClientRect(); return {t:e.innerText.slice(0,40), x:r.x,y:r.y,w:r.width,h:r.height, html:e.outerHTML.slice(0,600)}})""")
        (OUT / "selects.json").write_text(json.dumps(info, indent=1))
        sel = [i for i in info if i["t"].startswith("Collection") and i["w"] > 0][0]
        boxes["db_collection"] = {"x": sel["x"], "y": sel["y"], "width": sel["w"], "height": sel["h"]}
        cx, cy = sel["x"] + sel["w"] / 2, sel["y"] + sel["h"] - 20
        page.mouse.click(cx, cy); page.wait_for_timeout(1800)
        shot(page, "db_1_open")
        opts = page.evaluate("""() => [...document.querySelectorAll('[role=option], li')].filter(e=>e.offsetParent).map(e => {const r=e.getBoundingClientRect(); return {t:e.innerText.slice(0,40), x:r.x,y:r.y,w:r.width,h:r.height}})""")
        (OUT / "options.json").write_text(json.dumps(opts, indent=1))
        pol = [o for o in opts if o["t"].strip().startswith("Polymer")][0]
        boxes["db_opt_polymer"] = {"x": pol["x"], "y": pol["y"], "width": pol["w"], "height": pol["h"]}
        page.mouse.click(pol["x"] + pol["w"] / 2, pol["y"] + pol["h"] / 2); page.wait_for_timeout(9000)
        shot(page, "db_2_polymer")
        bb(page, "db_count2", page.get_by_text("matching records"))
        bb(page, "db_table2", page.get_by_test_id("stDataFrame"))
        search.click(); search.type("PCDTBT", delay=110); search.press("Enter"); page.wait_for_timeout(9000)
        shot(page, "db_3_query")
        bb(page, "db_count3", page.get_by_text("matching record"))
        bb(page, "db_table3", page.get_by_test_id("stDataFrame"))
    step("database", database)

    def record(page, rid, name, viewer):
        go(page, f"?page=database&record={rid}"); shot(page, name + "_0")
        for k, t, ex in [("exp", "Exp. Eg (eV)", True), ("homo", "HOMO (eV)", True), ("eg", "Eg (eV)", True), ("lot", "Level of theory", False),
                         ("prov", "Provenance", True), ("id", rid, True), ("ident", "Identity", True), ("props", "Electronic properties", True)]:
            bb(page, f"{name}_{k}", page.get_by_text(t, exact=ex))
        if viewer:
            tg = page.get_by_text("Show interactive 3D structure")
            bb(page, f"{name}_toggle", tg)
            tg.first.click(); page.wait_for_timeout(12000)
            page.evaluate("window.scrollTo(0,0)"); page.wait_for_timeout(500)
            shot(page, name + "_1_3d")
            bb(page, f"{name}_viewer", page.locator("iframe").last)
    step("pol", lambda pg: record(pg, "POL_0079", "pol", True))
    step("lit", lambda pg: record(pg, "LIT_0214", "lit", False))

    def stats(page):
        go(page, "?page=statistics"); shot(page, "stats")
        for k, t in [("snapshot", "Release snapshot"), ("dist", "Distribution overview")]:
            bb(page, "st_" + k, page.get_by_text(t))
    step("stats", stats)

    def ketcher(page):
        go(page, "?page=database")
        exp = page.locator("[data-testid='stExpander']").filter(has_text="Structure search").locator("summary")
        bb(page, "k_expander", exp)
        exp.first.click(); page.wait_for_timeout(2500)
        shot(page, "k_0")
        tg = page.get_by_text("Draw structure", exact=False)
        bb(page, "k_toggle", tg)
        tg.first.click(); page.wait_for_timeout(15000)
        shot(page, "ketcher")
        bb(page, "k_frame", page.locator("iframe").first)
    step("ketcher", ketcher)
    br.close()
(OUT / "boxes.json").write_text(json.dumps(boxes, indent=1))
(OUT / "log.txt").write_text("\n".join(log))
print("\n".join(log))
