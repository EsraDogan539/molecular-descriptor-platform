from pathlib import Path
from playwright.sync_api import sync_playwright
out = Path("docs/demo"); out.mkdir(parents=True, exist_ok=True)
B = "https://chalmoldb.streamlit.app/~/+"
log = []
def wake(page):
    page.goto(B + "/?page=home", timeout=120000)
    w = page.get_by_text("Yes, get this app back up")
    if w.count():
        w.first.click()
    page.wait_for_timeout(60000)
def shot(page, name, full=True):
    page.screenshot(path=out / f"{name}.png", full_page=False)
    if full:
        # streamlit scrolls inside a container; capture tall viewport separately
        pass
    log.append(name)
with sync_playwright() as p:
    b = p.chromium.launch()
    # tall captures (whole page) at 1600 wide, 1.5x
    ctx = b.new_context(viewport={"width": 1600, "height": 900}, device_scale_factor=1.5)
    page = ctx.new_page()
    wake(page)
    for name, q in [("home", "?page=home"), ("database", "?page=database"), ("statistics", "?page=statistics"),
                    ("analyze", "?page=analyze"), ("about", "?page=about"), ("documentation", "?page=documentation")]:
        page.goto(B + "/" + q, timeout=120000); page.wait_for_timeout(25000)
        shot(page, name + "_top")
        page.set_viewport_size({"width": 1600, "height": 4200}); page.wait_for_timeout(6000)
        page.screenshot(path=out / f"{name}_tall.png")
        page.set_viewport_size({"width": 1600, "height": 900}); page.wait_for_timeout(2000)
    # record pages with 3D
    for rid in ["EXT_0218", "EXT_0183", "DEV_0481"]:
        page.set_viewport_size({"width": 1600, "height": 3000})
        page.goto(B + f"/?page=database&record={rid}", timeout=120000); page.wait_for_timeout(25000)
        t = page.get_by_text("Show interactive 3D structure")
        if t.count():
            t.first.click(); page.wait_for_timeout(12000)
        page.screenshot(path=out / f"record_{rid}_tall.png")
        (out / f"record_{rid}.txt").write_text(page.locator("body").inner_text())
    # structure search with Ketcher
    page.set_viewport_size({"width": 1600, "height": 3000})
    page.goto(B + "/?page=database", timeout=120000); page.wait_for_timeout(25000)
    try:
        page.get_by_text("Structure search").first.click(); page.wait_for_timeout(3000)
        page.get_by_text("Draw structure").first.click(); page.wait_for_timeout(25000)
    except Exception as e:
        log.append(f"ketcher err {e}")
    page.screenshot(path=out / "ketcher_tall.png")
    ctx.close()
    # video of 3D viewer interaction
    ctx = b.new_context(viewport={"width": 1600, "height": 900}, device_scale_factor=1,
                        record_video_dir=str(out / "video"), record_video_size={"width": 1600, "height": 900})
    page = ctx.new_page()
    page.goto(B + "/?page=database&record=EXT_0218", timeout=120000); page.wait_for_timeout(30000)
    t = page.get_by_text("Show interactive 3D structure")
    if t.count():
        t.first.click(); page.wait_for_timeout(12000)
        fr = None
        for f in page.frames:
            try:
                if f.locator("canvas").count():
                    fr = f
            except Exception:
                pass
        if fr is not None:
            c = fr.locator("canvas").first
            c.scroll_into_view_if_needed(); page.wait_for_timeout(2000)
            bb = c.bounding_box(); log.append(f"canvas {bb}")
            cx, cy = bb["x"] + bb["width"] / 2, bb["y"] + bb["height"] / 2
            page.mouse.move(cx, cy); page.mouse.down()
            for i in range(160):
                page.mouse.move(cx + 220 * __import__('math').sin(i / 25), cy + 40 * __import__('math').sin(i / 40)); page.wait_for_timeout(40)
            page.mouse.up(); page.wait_for_timeout(1500)
    ctx.close()
    b.close()
(out / "log.txt").write_text("\n".join(map(str, log)) + "\n")
