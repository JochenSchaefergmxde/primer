import re
with open(r"C:\Users\internet\.gemini\antigravity\brain\e256fae4-44d0-4f26-a85b-f578875bdea1\.system_generated\steps\20\content.md", "r", encoding="utf-8") as f:
    data = f.read()

# Find strings that might be URLs
for m in set(re.findall(r'\"([a-zA-Z0-9_\-\./]+)\"', data)):
    if 'organism' in m.lower() or 'download' in m.lower():
        print("MATCH:", m)
