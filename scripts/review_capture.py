from pathlib import Path
from playwright.sync_api import sync_playwright
out = Path("docs/review"); out.mkdir(parents=True, exist_ok=True)
B = "https://chalmoldb.streamlit.app/~/+"
pages = [("home","?page=home"),("record","?page=database&record=DEV_0001"),("polymer","?page=database&record=EXT_0130"),
         ("statistics","?page=statistics"),("about","?page=about")]
log=[]
with sync_playwright() as p:
    b = p.chromium.launch(); page = b.new_page(viewport={"width": 1440, "height": 6000})
    page.goto(B+"/?page=home", timeout=120000)
    w = page.get_by_text("Yes, get this app back up")
    if w.count(): w.first.click()
    page.wait_for_timeout(240000)
    for name, q in pages:
        page.goto(B+"/"+q, timeout=120000); page.wait_for_timeout(25000)
        if name == "record":
            t = page.get_by_text("Show interactive 3D structure")
            if t.count(): t.first.click(); page.wait_for_timeout(8000)
        page.screenshot(path=out / f"{name}.png")
        txt = page.locator("body").inner_text()
        (out / f"{name}.txt").write_text(txt)
        log.append(f"{name}: {len(txt)} chars; exception={'Traceback' in txt or 'Error' in txt}")
    b.close()
(out / "log.txt").write_text("\n".join(log) + "\n")
