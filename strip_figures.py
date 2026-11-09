#!/usr/bin/env python3
"""notebook から画像の出力だけを落とした写しを作る。

6_pathway.ipynb は 95 MB あり、そのうち 97 MB 相当（base64）が画像なので git では
管理していない。代わりに、この写し（図なし版）を管理して GitHub でも中身を追えるようにする。

★ 落とすのは画像の出力だけ。cat() の文字出力・表（text/html）・エラーはすべて残すので、
  数値は読める。図は HTML（oligo_eptA_analysis.html）と /mnt/scratch/.../pathway/ にある。
★ 元の notebook には触らない。

使い方: python3 strip_figures.py <入力.ipynb> <出力.ipynb>
"""
import json, sys
from pathlib import Path

DROP = ("image/png", "image/jpeg", "image/svg+xml", "application/pdf")

NOTE = """> ★ **これは図を落とした写しです。**元の `6_pathway.ipynb` は 95 MB（うち画像 227 枚・97 MB 相当）あり、
> git では管理していません。`strip_figures.py` が自動で作ります（`.githooks/pre-commit`）。
>
> - **文字出力・表・エラーはすべて残しています。**数値はこのまま読めます
> - **図を見るには** `oligo_eptA_analysis.html`（まとめ HTML）か、
>   `/mnt/scratch/tanoue/261016_Xiaoxi_eptAoligo_2wFeces_Metatranscriptome/pathway/` を開いてください
> - **この写しを編集しないでください。**元の notebook から作り直されます
"""


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    src, dst = Path(sys.argv[1]), Path(sys.argv[2])
    nb = json.loads(src.read_text())

    n_img = 0
    for c in nb["cells"]:
        if c.get("cell_type") != "code":
            continue
        keep = []
        for o in c.get("outputs", []):
            d = o.get("data")
            if d:
                for k in list(d):
                    if k in DROP:
                        del d[k]; n_img += 1
                # 画像だけの出力は、残った text/plain が '<Figure ...>' のような
                # 中身の無いものになるので、丸ごと落とす
                if not d or set(d) <= {"text/plain"} and all(
                        "".join(v if isinstance(v, list) else [v]).strip().startswith(("<Figure", "<ggplot"))
                        for v in d.values()):
                    continue
            keep.append(o)
        c["outputs"] = keep

    # 先頭に断り書きを置く（元の見出しセルの前）
    nb["cells"].insert(0, {"cell_type": "markdown", "id": "nofig-note",
                           "metadata": {}, "source": NOTE.splitlines(keepends=True)})
    dst.write_text(json.dumps(nb, ensure_ascii=False, indent=1) + "\n")
    print(f"  {src.name} {src.stat().st_size/1e6:.1f} MB -> {dst.name} "
          f"{dst.stat().st_size/1e6:.1f} MB（画像の出力 {n_img} 個を削除）")


if __name__ == "__main__":
    main()
