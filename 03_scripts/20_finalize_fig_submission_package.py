#!/usr/bin/env python3
"""Create a clean upload-only directory after DOCX rendering and verification."""
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "06_manuscript" / "FIG_submission"
DST = SRC / "UPLOAD_PACKAGE"
if DST.exists():
    shutil.rmtree(DST)
(DST / "Figures").mkdir(parents=True, exist_ok=True)

files = {
    "01_Main_Manuscript_FIG.docx": "Main_Manuscript.docx",
    "02_Title_Page_FIG.docx": "Title_Page.docx",
    "03_Cover_Letter_FIG.docx": "Cover_Letter.docx",
    "04_Supplementary_Information_FIG.pdf": "ESM_1_Supplementary_Information.pdf",
    "06_Supplementary_Tables_FIG.xlsx": "ESM_2_Supplementary_Tables.xlsx",
}
for source, target in files.items():
    p = SRC / source
    if not p.exists():
        raise FileNotFoundError(p)
    shutil.copy2(p, DST / target)

for p in sorted((SRC / "Figures").glob("*.tif")):
    shutil.copy2(p, DST / "Figures" / p.name)

archive = SRC / "Functional_Integrative_Genomics_UPLOAD_PACKAGE"
shutil.make_archive(str(archive), "zip", root_dir=DST)
print(f"Upload package: {DST}")
print(f"Archive: {archive}.zip")
