from pathlib import Path
from playwright.sync_api import sync_playwright
out = Path("docs/live_check"); out.mkdir(parents=True, exist_ok=True)
lines = []
with sync_playwright() as p:
    b = p.chromium.launch(); page = b.new_page(viewport={"width": 1440, "height": 3000})
    for pg in ("home", "about", "statistics"):
        page.goto(f"https://chalmoldb.streamlit.app/~/+/?page={pg}", timeout=120000)
        wake = page.get_by_text("Yes, get this app back up")
        if wake.count(): wake.first.click()
        page.wait_for_timeout(40000)
        body = page.locator("body").inner_text()
        lines.append(f"{pg}: 3,488={'3,488' in body} 3,248={'3,248' in body} Developed-by={'Developed by Esra' in body} Team={'Hakan Kay' in body}")
    b.close()
(out / "numbers.txt").write_text("\n".join(lines) + "\n")
