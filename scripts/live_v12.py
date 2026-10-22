from pathlib import Path
from playwright.sync_api import sync_playwright
out = Path("docs/live_check"); out.mkdir(parents=True, exist_ok=True)
lines = []
with sync_playwright() as p:
    b = p.chromium.launch(); page = b.new_page(viewport={"width": 1440, "height": 2400})
    for rid in ("POL_0001", "POL_1200", "LIT_0001", "DEV_0001", "EXT_0300"):
        page.goto(f"https://chalmoldb.streamlit.app/~/+/?page=database&record={rid}", timeout=180000)
        wake = page.get_by_text("Yes, get this app back up")
        if wake.count(): wake.first.click()
        page.wait_for_timeout(45000)
        body = page.locator("body").inner_text()
        keys = [rid, "Experimental", "TD-DFT", "CCSD", "Nanoscale", "Kumar", "Traceback", "Error"]
        lines.append(rid + ": " + " ".join(f"{k}={k in body}" for k in keys))
        page.screenshot(path=str(out / f"rec_{rid}.png"), full_page=True)
        (out / f"rec_{rid}.txt").write_text(body)
    b.close()
(out / "v12_records.txt").write_text("\n".join(lines) + "\n")
