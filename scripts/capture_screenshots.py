"""Capture web-platform screenshots for the manuscript (Figure 4) from a locally running app.

    streamlit run app_v08.py --server.headless true --server.port 8501 &
    python scripts/capture_screenshots.py --url http://localhost:8501 --out docs/screenshots
"""

import argparse
from pathlib import Path

from playwright.sync_api import sync_playwright

RECORD = "EXT_0162"  # S/Te/S chalcogendiazoloquinoxaline monomer: exact structure + 3D coordinates


def wait_ready(page, text, timeout=120_000):
    page.get_by_text(text, exact=False).first.wait_for(timeout=timeout)
    page.wait_for_timeout(2500)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--url", default="http://localhost:8501")
    ap.add_argument("--out", default="docs/screenshots")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    log = []

    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page(viewport={"width": 1440, "height": 1000}, device_scale_factor=2)

        page.goto(f"{args.url}/?page=home")
        wait_ready(page, "Database at a glance")
        page.screenshot(path=out / "01_home.png")

        page.goto(f"{args.url}/?page=database")
        wait_ready(page, "matching record")
        page.get_by_text("matching record").first.scroll_into_view_if_needed()
        page.wait_for_timeout(1000)
        page.screenshot(path=out / "02_database_search.png")

        page.goto(f"{args.url}/?page=database&record={RECORD}")
        wait_ready(page, "matching record")
        heading = page.get_by_text("Inspect a record").first
        heading.scroll_into_view_if_needed()
        page.wait_for_timeout(1500)
        page.screenshot(path=out / "03_record_detail.png")
        try:
            page.get_by_text("Show interactive 3D structure").first.click()
            page.wait_for_timeout(8000)
            page.get_by_text("Drag to rotate").first.scroll_into_view_if_needed()
            page.wait_for_timeout(1500)
            page.screenshot(path=out / "04_record_3d.png")
        except Exception as exc:  # keep the other screenshots if the viewer cannot be captured
            log.append(f"3D viewer: {exc}")

        page.goto(f"{args.url}/?page=database")
        wait_ready(page, "matching record")
        try:
            page.get_by_text("Structure search").first.click()
            page.wait_for_timeout(1500)
            page.get_by_text("Draw structure").first.click()
            page.wait_for_timeout(15000)
            page.get_by_text("Structure search").first.scroll_into_view_if_needed()
            page.wait_for_timeout(1000)
            page.screenshot(path=out / "05_structure_editor.png")
        except Exception as exc:
            log.append(f"editor: {exc}")

        page.goto(f"{args.url}/?page=statistics")
        page.wait_for_timeout(12000)
        page.screenshot(path=out / "06_statistics.png")

        page.goto(f"{args.url}/?page=documentation")
        wait_ready(page, "Data access")
        page.get_by_text("Data access and programmatic use").first.scroll_into_view_if_needed()
        page.wait_for_timeout(1000)
        page.screenshot(path=out / "07_data_access.png")
        browser.close()

    (out / "capture_log.txt").write_text("\n".join(log) + "\n")
    print("done", log)


if __name__ == "__main__":
    main()
