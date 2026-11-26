import json
import csv

with open(r"C:\Users\internet\Downloads\repository\primer\organisms.json", "r", encoding="utf-8") as f:
    data = json.load(f)

csv_path = r"C:\Users\internet\Downloads\best_primers_db.csv"
with open(csv_path, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["Category", "Group", "OrgID", "OrgName", "OrgVersion", "OrgKey", "BestPrimersAvailable"])
    
    for category, groups in data.items():
        if isinstance(groups, dict):
            for group_name, items in groups.items():
                if isinstance(items, list):
                    for item in items:
                        if item.get("bestPrimers") == 1:
                            writer.writerow([
                                category,
                                group_name,
                                item.get("orgID", ""),
                                item.get("orgName", ""),
                                item.get("orgVersion", ""),
                                item.get("orgKey", ""),
                                item.get("bestPrimers", 0)
                            ])

print(f"Created CSV at {csv_path}")
