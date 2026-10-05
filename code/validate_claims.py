import json
from pathlib import Path

root = Path(__file__).resolve().parent.parent
data = json.loads((root / "code" / "claims.json").read_text(encoding="utf-8"))
claims = data["claims"]
kinds = {"paper_value", "paper_figure", "archived_replay", "illustrative_example"}
assert isinstance(claims, list) and claims
assert len({item["id"] for item in claims}) == len(claims)
for item in claims:
    assert item["kind"] in kinds, item["id"]
    source = root / item["source"]
    assert source.is_file(), item["source"]
print(f"OK: {len(claims)} claims; sources exist")
