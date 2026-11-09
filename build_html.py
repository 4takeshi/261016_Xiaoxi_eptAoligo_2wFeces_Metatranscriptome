#!/usr/bin/env python3
"""ステージ0〜7 の notebook 8本を1つの HTML にまとめる。

★ 2種類を必ず同時に作る（どちらかだけ古くなる事故を防ぐため）
    oligo_eptA_analysis.html        縮小版。GitHub で管理する
    oligo_eptA_analysis_full.html   原寸版。ローカルのみ（.gitignore）

手順
 1. 各 notebook を、保存済みの出力のまま quarto で markdown に変換する（再実行はしない）
 2. 見出しを1段下げ、「# 3. ...」などの章の下に入れる
 3. 1つの .qmd にまとめ、quarto で HTML にする（embed-resources で画像を埋め込む）

★ 縮小版は、.qmd の中の base64 画像を差し替えてから render する。
  notebook 本体には手を入れない。
  6_pathway.ipynb だけで画像 227 枚・72.7 MB あり、原寸だと HTML が
  115〜120 MB になって GitHub の 1ファイル 100 MB 上限を超える（実測）。

使い方: python3 build_html.py
"""
import base64, io, os, re, shutil, subprocess, sys
from pathlib import Path

project_dir = Path(__file__).resolve().parent
build_dir   = project_dir / "_build"
quarto      = Path.home() / "quarto" / "bin" / "quarto"

OUT_SMALL = project_dir / "oligo_eptA_analysis.html"
OUT_FULL  = project_dir / "oligo_eptA_analysis_full.html"

MAXW    = 1600   # 縮小版での画像の横幅の上限（px）
JPEG_Q  = 70     # 縮小版の JPEG 品質
#   2026-11-09: 88 -> 70 に下げた。セル3-4 が図1-3 の 38 枚を notebook に
#   出すようになり、88 のままだと縮小版が約 97 MB（GitHub 上限 100 MB の直前）
#   になるため。横幅 1600 px は維持したので、ブラウザで拡大したときの
#   遺伝子ラベルの読みやすさは変わらない。原寸版（_full）は無圧縮のまま

# (元ファイル, 章タイトル)
parts = [
    ("1_genome_copy.ipynb",            "ステージ0: ゲノムの収集"),
    ("2_reference_qc.ipynb",           "ステージ1: リファレンス QC とゲート1"),
    ("3_build_reference.ipynb",        "ステージ2: リファレンス構築"),
    ("4_reads.ipynb",                  "ステージ3-4: サンプル処理・QC・ゲート2"),
    ("5_1_DESeq2_A_community.ipynb",   "ステージ5 (A): 全株プール DESeq2"),
    ("5_2_DESeq2_B_perstrain.ipynb",   "ステージ5 (B): 株ごと独立 DESeq2"),
    ("6_pathway.ipynb",                "ステージ6: 機能解析"),
    ("7_genome_position.ipynb",        "ステージ7: ゲノム座位からの解析"),
]

front_matter = """---
title: 'eptA oligo-colonization gnotobiotic mouse metatranscriptome'
subtitle: '4株 defined community（Bth_VPI / Bth_St8 / Bun_8492 / Odo_St23）＋ pME008'
author: Takeshi Tanoue
date: last-modified
lang: ja
format:
  html:
    toc: true
    toc-depth: 3
    code-fold: true
    embed-resources: true
    css: styles.css
---

この HTML は、下の notebook を `build_html.py` で1つにまとめたものである。
**各 notebook は保存済みの出力をそのまま使っており、HTML 作成時に再実行はしていない。**

| 章 | 元ファイル |
|---|---|
""" + "".join(f"| {i}. {t} | `{s}` |\n" for i, (s, t) in enumerate(parts, start=1)) + """

解析の設計・決定事項は `CLAUDE.md` に、HTML 化と git の手順は `SETUP.md` にある。
"""


def copy_notebook(src, dst):
    """notebook を _build/ に写す。そのとき markdown セルの水平線 '---' を '***' に替える。
    ★ quarto は markdown 中の '---' ... '---' を YAML メタデータブロックとして読もうとし、
      'YAMLException: end of the stream or a document separator is expected' で落ちる（実測）。
      '***' も markdown の水平線なので見た目は変わらない。元の notebook には触らない。"""
    import json
    nb = json.loads(Path(src).read_text())
    n = 0
    for c in nb["cells"]:
        if c.get("cell_type") != "markdown":
            continue
        out = []
        for line in "".join(c["source"]).split("\n"):
            if line.strip() in ("---", "___"):
                line = "***"; n += 1
            out.append(line)
        o = "\n".join(out)
        c["source"] = [l + "\n" for l in o.split("\n")[:-1]] + [o.split("\n")[-1]]
    Path(dst).write_text(json.dumps(nb, ensure_ascii=False))
    if n:
        print(f"  水平線 '---' を '***' に置換: {n} 本")


def notebook_to_markdown(src):
    """notebook を quarto で markdown に変換し、見出しを1段下げた本文を返す。"""
    md_name = Path(src).with_suffix(".md").name
    subprocess.run([str(quarto), "render", src, "--to", "markdown", "--output", md_name],
                   cwd=build_dir, check=True)
    text = (build_dir / md_name).read_text()

    # quarto は最初の見出しを title として YAML に移すので、本文の先頭に戻す
    body, title = text, None
    m = re.match(r"---\n(.*?)\n---\n", text, re.S)
    if m:
        body = text[m.end():]
        t = re.search(r"^title:\s*(.*)$", m.group(1), re.M)
        if t:
            title = t.group(1).strip().strip("'\"")

    # ★ コードブロックの外にある見出しだけを1段下げる
    lines, in_fence = [], False
    for line in body.splitlines():
        if line.lstrip().startswith("```"):
            in_fence = not in_fence
        elif not in_fence and re.match(r"^#{1,5} ", line):
            line = "#" + line
        lines.append(line)
    return (f"## {title}\n\n" if title else "") + "\n".join(lines) + "\n"


def shrink_images(qmd_text):
    """縮小版のために、_build/*_files/ の PNG を縮小した JPEG に作り直し、
    .qmd の参照を .png -> .jpg に書き換えた本文を返す。

    ★ quarto は --to markdown のとき、画像を base64 ではなく
      <notebook>_files/figure-markdown/*.png というファイルに書き出す。
      そのため .qmd の文字列を置換しても画像は小さくならない（実測で 0 枚）。
      ファイルのほうを作り直す。
    ★ 元の notebook にも、原寸版の HTML にも影響しない
      （原寸版はこの関数を呼ぶ前に render 済み）。"""
    from PIL import Image
    n = a = b = 0
    for png in sorted(build_dir.rglob("*.png")):
        a += png.stat().st_size
        try:
            im = Image.open(png); im.load()
        except Exception:
            continue
        if im.width > MAXW:
            h = max(1, round(im.height * MAXW / im.width))
            im = im.resize((MAXW, h), Image.LANCZOS)
        if im.mode in ("RGBA", "LA", "P"):          # JPEG は透過を持てない
            bg = Image.new("RGB", im.size, "white")
            im = im.convert("RGBA"); bg.paste(im, mask=im.split()[-1]); im = bg
        else:
            im = im.convert("RGB")
        jpg = png.with_suffix(".jpg")
        im.save(jpg, "JPEG", quality=JPEG_Q, optimize=True)
        b += jpg.stat().st_size; n += 1
    print(f"  画像を縮小: {n} 枚 / {a/1e6:.1f} MB -> {b/1e6:.1f} MB")
    return re.sub(r"(\([^)]*_files/[^)]*)\.png\)", r"\1.jpg)", qmd_text)


# ★ GitHub は 1ファイル 100 MB で push を拒否する。縮小版がそこに近づいたら知らせる。
#   超えたら build_html.py の MAXW（横幅の上限）と JPEG_Q（品質）を下げる。
#   実測の目安: 1600px/q88 -> 79 MB / 1400px/q85 -> 57 MB / 1200px/q82 -> 45 MB
GH_LIMIT_MB = 100
GH_WARN_MB  = 90


def render(qmd_name, out_path, label, check_limit=False):
    subprocess.run([str(quarto), "render", qmd_name, "--to", "html"], cwd=build_dir, check=True)
    built = build_dir / Path(qmd_name).with_suffix(".html").name
    shutil.copy(built, out_path)
    mb = os.path.getsize(out_path) / 1e6
    print(f"■ {label}: {out_path.name}  {mb:.1f} MB")
    if check_limit:
        if mb >= GH_LIMIT_MB:
            print(f"!! {mb:.1f} MB は GitHub の 1ファイル {GH_LIMIT_MB} MB を超えています。"
                  f"push が拒否されます。MAXW / JPEG_Q を下げて作り直してください")
        elif mb >= GH_WARN_MB:
            print(f"!! {mb:.1f} MB。GitHub の上限 {GH_LIMIT_MB} MB まで残り {GH_LIMIT_MB-mb:.1f} MB です")
        elif mb >= 50:
            print(f"   （GitHub は 50 MB 超で警告を出します。更新のたびに約 {mb:.0f} MB "
                  f"リポジトリに積まれる点に注意）")


def main():
    if not quarto.exists():
        sys.exit(f"quarto が見つかりません: {quarto}")
    if build_dir.exists():
        shutil.rmtree(build_dir)
    build_dir.mkdir()

    sections = [front_matter]
    for i, (src, title) in enumerate(parts, start=1):
        p = project_dir / src
        if not p.exists():
            print(f"!! {src} がありません。飛ばします"); continue
        print(f"[{i}/{len(parts)}] {src}")
        copy_notebook(p, build_dir / src)
        sections.append(f"# {i}. {title}\n\n元ファイル: `{src}`\n")
        sections.append(notebook_to_markdown(src))

    shutil.copy(project_dir / "styles.css", build_dir / "styles.css")
    qmd_full = build_dir / "oligo_eptA_analysis_full.qmd"
    qmd_full.write_text("\n".join(sections))
    print(f"  .qmd {qmd_full.stat().st_size/1e6:.1f} MB")

    # ★ 必ず2本とも作る。片方だけ古くなる事故を防ぐ
    render(qmd_full.name, OUT_FULL, "原寸版（ローカルのみ）")
    qmd_small = build_dir / "oligo_eptA_analysis.qmd"
    qmd_small.write_text(shrink_images(qmd_full.read_text()))
    render(qmd_small.name, OUT_SMALL, "縮小版（GitHub 管理）", check_limit=True)


if __name__ == "__main__":
    main()
