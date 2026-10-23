"""Open the public ChalMolDB app, wake it if Streamlit put it to sleep, and check that it renders.

Exit code 1 (failed workflow run, e-mailed to the repository owner) when the app shows an
error page or the home page metrics do not appear.
"""
import sys
from playwright.sync_api import sync_playwright

URL = "https://chalmoldb.streamlit.app/"
EXPECTED = ["ChalMolDB", "Records"]
ERRORS = ["supposed to happen", "Traceback", "TypeError", "Error running app"]


def check(page):
    page.goto(URL, timeout=180_000)
    page.wait_for_timeout(5_000)
    wake = page.get_by_text("Yes, get this app back up")
    if wake.count():
        print("App was asleep; waking it")
        wake.first.click()
        page.wait_for_timeout(90_000)
    page.wait_for_timeout(30_000)
    frame_text = page.locator("body").inner_text()
    for f in page.frames:
        try:
            frame_text += "\n" + f.locator("body").inner_text(timeout=5_000)
        except Exception:
            pass
    bad = [e for e in ERRORS if e in frame_text]
    missing = [e for e in EXPECTED if e not in frame_text]
    return bad, missing


with sync_playwright() as p:
    browser = p.chromium.launch()
    page = browser.new_page(viewport={"width": 1280, "height": 900})
    for attempt in (1, 2):
        bad, missing = check(page)
        if not bad and not missing:
            print("OK: app is up")
            browser.close()
            sys.exit(0)
        print(f"Attempt {attempt}: errors={bad} missing={missing}")
    browser.close()
    sys.exit(1)
