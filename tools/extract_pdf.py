"""Extract the text layer of the source PDF into per-page files and a combined file.

This only reads publicly rendered text (definitions, equations as typeset,
numerical tables). It does not touch any program source code.
"""
import os
import fitz  # PyMuPDF

PDF = r"source\Questioning Fundamental Principles of Organic Chemistry.PDF"
OUT_DIR = os.path.join("data", "book_text")
os.makedirs(OUT_DIR, exist_ok=True)

doc = fitz.open(PDF)
n = doc.page_count

combined_lines = []
stats = []
for pno in range(n):
    page = doc[pno]
    try:
        txt = page.get_text("text")
    except Exception as e:  # noqa: BLE001
        txt = ""
        print(f"[warn] page {pno} text extraction failed: {e}")
    try:
        imgs = page.get_images(full=True)
    except Exception:
        imgs = []
    stats.append((pno, len(txt), len(imgs)))
    # per-page file (1-indexed for human readability)
    with open(os.path.join(OUT_DIR, f"page_{pno+1:03d}.txt"), "w", encoding="utf-8") as f:
        f.write(txt)
    combined_lines.append(f"\n\n===== PAGE {pno+1} (chars={len(txt)}, images={len(imgs)}) =====\n")
    combined_lines.append(txt)

with open(os.path.join("data", "book_full.txt"), "w", encoding="utf-8") as f:
    f.writelines(combined_lines)

total_chars = sum(s[1] for s in stats)
total_imgs = sum(s[2] for s in stats)
print(f"pages={n} total_chars={total_chars} total_images={total_imgs}")
low = [s for s in stats if s[1] < 50]
print(f"pages_with_little_text(<50 chars)={len(low)}: {[s[0]+1 for s in low][:40]}")
