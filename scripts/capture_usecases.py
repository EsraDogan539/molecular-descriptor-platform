"""Temporary: preview screenshots of the use-case section, quick start and a preset result."""
from pathlib import Path
from playwright.sync_api import sync_playwright
out = Path("docs/preview"); out.mkdir(parents=True, exist_ok=True)
with sync_playwright() as p:
    b = p.chromium.launch(); page = b.new_page(viewport={"width": 1440, "height": 1000})
    page.goto("http://localhost:8501/?page=home"); page.get_by_text("What you can do with ChalMolDB").first.wait_for(timeout=120000)
    page.wait_for_timeout(2500); page.get_by_text("What you can do with ChalMolDB").first.scroll_into_view_if_needed()
    page.wait_for_timeout(800); page.screenshot(path=out / "home_usecases.png")
    page.goto("http://localhost:8501/?page=documentation"); page.get_by_text("Quick start").first.wait_for(timeout=120000)
    page.wait_for_timeout(2000); page.get_by_text("Quick start").first.scroll_into_view_if_needed(); page.wait_for_timeout(800)
    page.screenshot(path=out / "quick_start.png")
    page.goto("http://localhost:8501/?page=database&preset=te_low_gap"); page.get_by_text("matching record").first.wait_for(timeout=120000)
    page.wait_for_timeout(2500); page.screenshot(path=out / "preset_te.png", full_page=False)
    b.close()
