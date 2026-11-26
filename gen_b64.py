import base64
with open(r"C:\Users\internet\Downloads\repository\primer\primer_suche.py", "rb") as f:
    b64 = base64.b64encode(f.read()).decode("ascii")

cmd = f"echo '{b64}' | base64 -d > /mnt/sdc1/scripts/primer_suche.py && chmod +x /mnt/sdc1/scripts/primer_suche.py"
with open(r"C:\Users\internet\Downloads\repository\primer\b64_cmd.txt", "w", encoding="utf-8") as f:
    f.write(cmd)
print("Saved b64_cmd.txt, length:", len(cmd))
