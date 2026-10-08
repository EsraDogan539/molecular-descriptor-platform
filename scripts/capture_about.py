"""Temporary: screenshots of the home header and the About page (team, citation)."""
from pathlib import Path
from playwright.sync_api import sync_playwright
out = Path("docs/screenshots_check"); out.mkdir(parents=True, exist_ok=True)
with sync_playwright() as p:
    b = p.chromium.launch()
    page = b.new_page(viewport={"width": 1440, "height": 1000}, device_scale_factor=1)
    page.goto("http://localhost:8501/?page=home"); page.get_by_text("Database at a glance").first.wait_for(timeout=120000)
    page.wait_for_timeout(3000); page.screenshot(path=out / "home.png")
    page.get_by_text("Developed by Esra Nur Doğan").first.scroll_into_view_if_needed(); page.wait_for_timeout(1000)
    page.screenshot(path=out / "home_footer.png")
    page.goto("http://localhost:8501/?page=about"); page.get_by_text("How to cite").first.wait_for(timeout=120000)
    page.wait_for_timeout(2000)
    page.get_by_text("Team").first.scroll_into_view_if_needed(); page.wait_for_timeout(800)
    page.screenshot(path=out / "about.png")
    b.close()
