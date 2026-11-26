import re

with open(r"C:\Users\internet\.gemini\antigravity\brain\e256fae4-44d0-4f26-a85b-f578875bdea1\.system_generated\steps\20\content.md", "r", encoding="utf-8") as f:
    data = f.read()

# Find URLs in the code
urls = re.findall(r'https?://[^\s\"\'\`]+', data)
print("HTTP URLs found:")
for u in set(urls):
    if 'api' in u or 'download' in u or 'biodb' in u:
        print(u)

paths = re.findall(r'\"(/[a-zA-Z0-9_-]+/[^\s\"\'\`]*)\"', data)
print("\nPaths found:")
for p in set(paths):
    if 'api' in p or 'download' in p:
        print(p)
