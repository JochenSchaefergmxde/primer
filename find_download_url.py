import re

with open(r"C:\Users\internet\.gemini\antigravity\brain\e256fae4-44d0-4f26-a85b-f578875bdea1\.system_generated\steps\359\content.md", "r", encoding="utf-8") as f:
    data = f.read()

print("Length:", len(data))

# Search for window.open or href or download or url constructions
matches = re.findall(r'.{0,100}(?:download|http|best|all).{0,100}', data, re.I)
for m in set(matches[:30]):
    print("MATCH:", m)
