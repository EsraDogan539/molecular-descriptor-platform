import json, urllib.request, time
from pathlib import Path
out = Path("docs/live_check"); out.mkdir(parents=True, exist_ok=True)
res = {}
for attempt in range(10):
    try:
        url = "https://zenodo.org/api/records?q=conceptrecid:22903549&all_versions=true&sort=mostrecent&size=5"
        d = json.load(urllib.request.urlopen(url, timeout=60))
        res = [{"id": h["id"], "doi": h.get("doi"), "title": h["metadata"]["title"], "version": h["metadata"].get("version"),
                "date": h["metadata"].get("publication_date"), "created": h.get("created")} for h in d["hits"]["hits"]]
        if any("1.2" in (r["title"] or "") for r in res): break
    except Exception as e:
        res = {"error": repr(e)}
    time.sleep(60)
(out / "zenodo_v12.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
print(json.dumps(res, indent=1, ensure_ascii=False))
