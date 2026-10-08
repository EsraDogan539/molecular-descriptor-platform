from pathlib import Path
from playwright.sync_api import sync_playwright
out = Path("docs/wp"); out.mkdir(parents=True, exist_ok=True)
B = "http://127.0.0.1:9400"
with sync_playwright() as p:
    b = p.chromium.launch(); pg = b.new_page(viewport={"width": 1280, "height": 1400})
    for name, url in [("home", "/"), ("menus", "/wp-admin/nav-menus.php"), ("pages", "/wp-admin/edit.php?post_type=page"), ("newpage", "/wp-admin/post-new.php?post_type=page")]:
        pg.goto(B + url, timeout=120000); pg.wait_for_timeout(8000)
        pg.screenshot(path=out / f"{name}.png")
        (out / f"{name}.txt").write_text(pg.locator("body").inner_text()[:4000])
    b.close()
