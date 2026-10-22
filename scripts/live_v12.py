from pathlib import Path
from playwright.sync_api import sync_playwright
out = Path("docs/live_check"); out.mkdir(parents=True, exist_ok=True)
lines = []
with sync_playwright() as p:
    b = p.chromium.launch(); page = b.new_page(viewport={"width": 1440, "height": 2400})
    for pg in ("home", "about", "database"):
        page.goto(f"https://chalmoldb.streamlit.app/~/+/?page={pg}", timeout=180000)
        wake = page.get_by_text("Yes, get this app back up")
        if wake.count(): wake.first.click()
        page.wait_for_timeout(60000)
        body = page.locator("body").inner_text()
        keys = ["5,029", "3,488", "4,229", "1,371", "Polymer", "Literature", "Nanoscale", "Database v1.2", "Database v1.1"]
        lines.append(pg + ": " + " ".join(f"{k}={k in body}" for k in keys))
        page.screenshot(path=str(out / f"v12_{pg}.png"), full_page=True)
        (out / f"v12_{pg}.txt").write_text(body)
    page.goto("https://chalmoldb.streamlit.app/~/+/?record=POL_0001", timeout=180000)
    page.wait_for_timeout(60000)
    body = page.locator("body").inner_text()
    lines.append("record POL_0001: " + str("POL_0001" in body) + " exp=" + str("Experimental" in body))
    page.screenshot(path=str(out / "v12_record.png"), full_page=True)
    b.close()
(out / "v12.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
