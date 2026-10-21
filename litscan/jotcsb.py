import urllib.request, re, json
from pathlib import Path
out = Path("docs/litscan/jotcsb"); out.mkdir(parents=True, exist_ok=True)
UA = {"User-Agent": "Mozilla/5.0 (ChalMolDB literature check; mailto:esra.6615@gmail.com)"}
log = {}
for name, url in [("article.pdf", "https://dergipark.org.tr/tr/download/article-file/286973"),
                  ("article_gov.pdf", "http://dergipark.gov.tr/download/article-file/286973")]:
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120) as r:
            data = r.read(); (out / name).write_bytes(data); log[name] = [r.geturl(), len(data), data[:8].decode("latin1")]
    except Exception as e:
        log[name] = str(e)
(out / "log.json").write_text(json.dumps(log, indent=1)); print(log)
