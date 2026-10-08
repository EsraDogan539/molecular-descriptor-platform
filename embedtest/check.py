import json, threading, http.server, functools
from pathlib import Path
from playwright.sync_api import sync_playwright
out = Path("docs/embed"); out.mkdir(parents=True, exist_ok=True)
H = functools.partial(http.server.SimpleHTTPRequestHandler, directory="embedtest")
srv = http.server.ThreadingHTTPServer(("127.0.0.1", 8765), H); threading.Thread(target=srv.serve_forever, daemon=True).start()
log = []
with sync_playwright() as p:
    b = p.chromium.launch()
    # direct visit headers
    ctx = b.new_context(); pg = ctx.new_page()
    def onresp(r):
        if "streamlit.app" in r.url or "share.streamlit.io" in r.url:
            h = r.headers
            log.append({"url": r.url[:140], "status": r.status, "xfo": h.get("x-frame-options"), "csp": (h.get("content-security-policy") or "")[:300], "loc": h.get("location")})
    pg.on("response", onresp)
    pg.goto("https://chalmoldb.streamlit.app/?embed=true&embedded=1", timeout=120000); pg.wait_for_timeout(90000)
    pg.screenshot(path=out / "direct.png")
    # embedded in a third-party page
    log.clear(); ctx2 = b.new_context(); pg2 = ctx2.new_page()
    def onresp2(r):
        if r.request.resource_type == "document":
            h = r.headers
            log.append({"url": r.url[:140], "status": r.status, "xfo": h.get("x-frame-options"), "csp": (h.get("content-security-policy") or "")[:300], "loc": h.get("location"), "setcookie": (h.get("set-cookie") or "")[:200]})
    pg2.on("response", onresp2)
    pg2.on("requestfailed", lambda q: log.append({"FAILED": q.url[:140], "err": q.failure}))
    msgs = []; pg2.on("console", lambda m: msgs.append(m.text[:200]))
    pg2.goto("http://127.0.0.1:8765/page.html", timeout=120000); pg2.wait_for_timeout(90000)
    pg2.screenshot(path=out / "embedded.png")
    frames = [{"url": f.url[:140]} for f in pg2.frames]
    (out / "log.json").write_text(json.dumps({"responses": log[:60], "frames": frames, "console": msgs[:40]}, indent=1))
    b.close()
