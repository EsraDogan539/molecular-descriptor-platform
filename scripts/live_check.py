"""Check the deployed app (temporary helper)."""
import sys, time
from pathlib import Path
from playwright.sync_api import sync_playwright

BASE = "https://chalmoldb.streamlit.app"
INNER = BASE + "/~/+"
out = Path("docs/live_check"); out.mkdir(parents=True, exist_ok=True)
log = []
def note(msg):
    print(msg); log.append(msg)

with sync_playwright() as p:
    b = p.chromium.launch()
    page = b.new_page(viewport={"width": 1440, "height": 1000})
    deadline = time.time() + 1500
    ok = False
    while time.time() < deadline:
        try:
            page.goto(INNER + "/?page=home", timeout=90_000)
            wake = page.get_by_text("Yes, get this app back up")
            if wake.count():
                wake.first.click(); note("woke app")
            page.get_by_text("Database at a glance").first.wait_for(timeout=180_000)
            page.wait_for_timeout(3000)
            if page.get_by_text("3,248").count():
                ok = True; break
            note("old version still served; waiting")
        except Exception as exc:
            note(f"home not ready: {type(exc).__name__}")
        time.sleep(60)
    note(f"new version live: {ok}")
    page.screenshot(path=out / "home.png")

    page.goto(INNER + "/?page=database&record=EXT_0162", timeout=90_000)
    page.get_by_text("Inspect a record").first.wait_for(timeout=180_000)
    page.wait_for_timeout(3000)
    note(f"deep link selects EXT_0162: {page.get_by_text('EXT_0162 — STeS').count() > 0}")
    try:
        page.get_by_text("Show interactive 3D structure").first.click()
        page.wait_for_timeout(10000)
        canvas = page.frame_locator("iframe").locator("canvas")
        note(f"3D viewer canvas present: {canvas.count() > 0}")
        page.get_by_text("Drag to rotate").first.scroll_into_view_if_needed()
        page.screenshot(path=out / "viewer.png")
    except Exception as exc:
        note(f"3D viewer failed: {exc}")
    note(f"SDF download button: {page.get_by_text('Download 3D structure (SDF)').count() > 0}")

    page.goto(INNER + "/?page=database", timeout=90_000)
    page.get_by_text("matching records").first.wait_for(timeout=180_000)
    try:
        page.get_by_text("Structure search").first.click(); page.wait_for_timeout(1500)
        page.get_by_text("Draw structure").first.click(); page.wait_for_timeout(15000)
        frames = [f.url for f in page.frames]
        note(f"Ketcher frame loaded: {any('ketcher' in u.lower() or 'component' in u.lower() for u in frames)}")
        page.screenshot(path=out / "editor.png")
    except Exception as exc:
        note(f"editor failed: {exc}")

    page.goto(INNER + "/?page=documentation", timeout=90_000)
    page.get_by_text("Data access and programmatic use").first.wait_for(timeout=180_000)
    note("documentation data-access section present: True")

    # Outer URL (as shared publicly) must forward the record parameter
    page.goto(BASE + "/?page=database&record=EXT_0162", timeout=90_000)
    page.wait_for_timeout(30000)
    found = any(f.locator("text=EXT_0162 — STeS").count() > 0 for f in page.frames)
    note(f"public permanent link works: {found}")
    page.screenshot(path=out / "public_link.png")
    b.close()
(out / "log.txt").write_text("\n".join(log) + "\n")
