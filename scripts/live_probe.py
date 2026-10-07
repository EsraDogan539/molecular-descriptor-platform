import time
from pathlib import Path
from playwright.sync_api import sync_playwright
out = Path("docs/live_check"); out.mkdir(parents=True, exist_ok=True)
lines = []
with sync_playwright() as p:
    b = p.chromium.launch()
    page = b.new_page(viewport={"width": 1440, "height": 1000},
                      user_agent="Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/130.0 Safari/537.36")
    for name, url in (("outer", "https://chalmoldb.streamlit.app/?page=home"), ("inner", "https://chalmoldb.streamlit.app/~/+/?page=home")):
        r = page.goto(url, timeout=120_000)
        lines.append(f"== {name} status {r.status if r else None} final {page.url}")
        for t in range(12):
            page.wait_for_timeout(10_000)
            for btn in ("Yes, get this app back up!", "Yes, get this app back up"):
                if page.get_by_text(btn).count():
                    page.get_by_text(btn).first.click(); lines.append("clicked wake")
            texts = []
            for f in page.frames:
                try:
                    texts.append((f.url[:90], f.locator("body").inner_text(timeout=5000)[:300].replace("\n", " | ")))
                except Exception as e:
                    texts.append((f.url[:90], f"ERR {type(e).__name__}"))
            if any("3,248" in t or "3,360" in t for _, t in texts):
                break
        for u, t in texts:
            lines.append(f"  frame {u}: {t}")
        page.screenshot(path=out / f"probe_{name}.png")
    b.close()
(out / "probe.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
