from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urldefrag, urlparse
ROOT = Path.cwd()
SKIP = ("http:", "https:", "//", "__SITE" + "_URL__", "mailto:", "data:", "javascript:")
class Scanner(HTMLParser):
    def __init__(self):
        super().__init__(); self.links = []
    def handle_starttag(self, _, attrs):
        self.links.extend(value for name, value in attrs if name in {"href", "src"} and value)
for html in ROOT.rglob("*.html"):
    parser = Scanner()
    parser.feed(html.read_text(encoding="utf-8"))
    for raw in parser.links:
        value, _ = urldefrag(raw)
        if not value or value.startswith(SKIP):
            continue
        current = html.parent
        for part in urlparse(value).path.replace("\\", "/").split("/"):
            if not part or part == ".":
                continue
            if part == "..":
                current = current.parent
                continue
            names = {entry.name for entry in current.iterdir()}
            assert part in names, f"case-sensitive missing target: {html} -> {raw}"
            current /= part
        assert current.is_file(), f"target is not a file: {html} -> {raw}"
print("OK: local HTML href/src targets exist with exact case")
