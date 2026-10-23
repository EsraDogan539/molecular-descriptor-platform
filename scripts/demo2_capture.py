"""Capture v1.2 site states for the demo video (2x DPR) plus element boxes."""
import json, traceback
from pathlib import Path
from playwright.sync_api import sync_playwright
OUT = Path("docs/demo2"); OUT.mkdir(parents=True, exist_ok=True)
BASE = "https://chalmoldb.streamlit.app/~/+/"
boxes, log = {}, []

def bb(page, key, loc):
    try:
        b = loc.first.bounding_box(timeout=8000)
        boxes[key] = b
    except Exception as e:
        boxes[key] = None; log.append(f"box {key}: {e!r}"[:300])

def shot(page, name, full=True):
    page.screenshot(path=str(OUT / f"{name}.png"), full_page=full)
    (OUT / f"{name}.txt").write_text(page.locator("body").inner_text())
    log.append(f"shot {name}")

def go(page, q, wait=45000):
    page.goto(BASE + q, timeout=180000)
    w = page.get_by_text("Yes, get this app back up")
    if w.count(): w.first.click()
    page.wait_for_timeout(wait)
    # hide streamlit chrome
    page.add_style_tag(content="header[data-testid='stHeader'],[data-testid='stToolbar'],[data-testid='stDecoration'],#MainMenu,footer,.stDeployButton,[data-testid='manage-app-button'],iframe[title='streamlitApp']+div{display:none!important}")
    page.wait_for_timeout(800)

def step(name, fn):
    try:
        fn(); log.append(f"ok {name}")
    except Exception:
        log.append(f"FAIL {name}: " + traceback.format_exc()[-800:])

with sync_playwright() as p:
    br = p.chromium.launch()
    ctx = br.new_context(viewport={"width": 1440, "height": 900}, device_scale_factor=2)
    page = ctx.new_page()

    def home():
        go(page, "?page=home")
        shot(page, "home")
        bb(page, "home_explore", page.get_by_text("Explore database", exact=True))
        bb(page, "home_glance", page.get_by_text("Database at a glance", exact=False))
        bb(page, "home_polymer_card", page.get_by_text("Polymer collection · 1,323 records"))
    step("home", home)

    def database():
        go(page, "?page=database")
        shot(page, "db_0")
        bb(page, "db_search", page.get_by_placeholder("ID, molecule/system, structure, method, source or reference"))
        sel = page.get_by_test_id("stSelectbox")
        bb(page, "db_collection", sel.nth(0))
        bb(page, "db_chalcogen", sel.nth(2))
        bb(page, "db_count", page.get_by_text("matching records"))
        bb(page, "db_table", page.get_by_test_id("stDataFrame"))
        sel.nth(0).click(); page.wait_for_timeout(1500)
        shot(page, "db_1_open", full=False)
        opt = page.get_by_role("option", name="Polymer")
        bb(page, "db_opt_polymer", opt)
        opt.first.click(); page.wait_for_timeout(9000)
        shot(page, "db_2_polymer")
        bb(page, "db_count2", page.get_by_text("matching records"))
        bb(page, "db_table2", page.get_by_test_id("stDataFrame"))
        inp = page.get_by_placeholder("ID, molecule/system, structure, method, source or reference")
        inp.click(); inp.type("PCDTBT", delay=120); inp.press("Enter"); page.wait_for_timeout(9000)
        shot(page, "db_3_query")
        bb(page, "db_count3", page.get_by_text("matching record"))
        bb(page, "db_table3", page.get_by_test_id("stDataFrame"))
    step("database", database)

    def record(rid, name, viewer=True):
        go(page, f"?page=database&record={rid}")
        shot(page, name + "_0")
        for k, t in [("exp", "Exp. Eg (eV)"), ("homo", "HOMO (eV)"), ("eg", "Eg (eV)"), ("lot", "Level of theory"), ("prov", "Provenance"), ("id", rid)]:
            bb(page, f"{name}_{k}", page.get_by_text(t, exact=(k != "lot")))
        if viewer:
            tg = page.get_by_text("Show interactive 3D structure")
            bb(page, f"{name}_toggle", tg)
            tg.first.click(); page.wait_for_timeout(12000)
            shot(page, name + "_1_3d")
            bb(page, f"{name}_viewer", page.locator("iframe").last)
    step("record_pol", lambda: record("POL_0079", "pol"))
    step("record_lit", lambda: record("LIT_0214", "lit", viewer=False))

    def stats():
        go(page, "?page=statistics")
        shot(page, "stats")
        for k, t in [("snapshot", "Release snapshot"), ("dist", "Distribution overview")]:
            bb(page, "st_" + k, page.get_by_text(t))
    step("stats", stats)

    def ketcher():
        go(page, "?page=database")
        page.get_by_text("Structure search").first.click(); page.wait_for_timeout(2500)
        tg = page.get_by_text("Draw structure")
        bb(page, "k_toggle", tg)
        tg.first.click(); page.wait_for_timeout(15000)
        shot(page, "ketcher")
        bb(page, "k_frame", page.locator("iframe").first)
    step("ketcher", ketcher)

    def docs():
        go(page, "?page=documentation")
        shot(page, "docs")
        bb(page, "docs_data", page.get_by_text("Data access", exact=False))
        go(page, "?page=about")
        shot(page, "about")
        bb(page, "about_cite", page.get_by_text("How to cite", exact=False))
    step("docs", docs)
    br.close()
(OUT / "boxes.json").write_text(json.dumps(boxes, indent=1))
(OUT / "log.txt").write_text("\n".join(log))
print("\n".join(log))
