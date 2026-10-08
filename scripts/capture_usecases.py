"""Temporary: preview the two hero texts."""
from pathlib import Path
from playwright.sync_api import sync_playwright
out = Path("docs/preview"); out.mkdir(parents=True, exist_ok=True)
with sync_playwright() as p:
    b = p.chromium.launch(); page = b.new_page(viewport={"width": 1440, "height": 1000}, device_scale_factor=1)
    for v in ("long", "short"):
        page.goto(f"http://localhost:8501/?page=home&hero={v}"); page.get_by_text("What you can do with ChalMolDB").first.wait_for(timeout=120000)
        page.wait_for_timeout(3000); page.screenshot(path=out / f"hero_{v}.png")
    b.close()
