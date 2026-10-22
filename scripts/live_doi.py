from pathlib import Path
from playwright.sync_api import sync_playwright
out = Path("docs/live_check"); out.mkdir(parents=True, exist_ok=True)
lines = []
with sync_playwright() as p:
    b = p.chromium.launch(); page = b.new_page(viewport={"width": 1440, "height": 2400})
    for pg in ("home", "about"):
        page.goto(f"https://chalmoldb.streamlit.app/~/+/?page={pg}", timeout=180000)
        wake = page.get_by_text("Yes, get this app back up")
        if wake.count(): wake.first.click()
        page.wait_for_timeout(50000)
        body = page.locator("body").inner_text()
        lines.append(pg + ": " + " ".join(f"{k}={k in body}" for k in ["23558298", "Database v1.2", "5,029"]))
    b.close()
(out / "doi.txt").write_text("\n".join(lines) + "\n")
