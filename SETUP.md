# HTML 化・git/GitHub 管理の構築記録（2026-11-09）

ステージ0〜7 がひと段落したため、`261013_Blastp_acdA` と同じように、解析の記録を
1つの HTML にまとめ、git/GitHub で管理する環境を Claude Code と一緒に作った。
このファイルは、その決定事項・行った作業・今後の保存手順をまとめたものである。

## 決定事項

**Q1. notebook 8本を1つの HTML にまとめられるか？**
→ **まとめられる。`build_html.py` で `oligo_eptA_analysis.html` を作る。**

- notebook が8本あり、Python と R が混ざっている。1つの notebook にはカーネルを1つしか
  指定できないので、notebook を直接つなげる方法は使えない。
- そこで各 notebook を**保存済みの出力のまま** quarto で markdown に変換し、1つの `.qmd` に
  まとめてから HTML にする。**再実行はしない**（bowtie2 や eggNOG-mapper が走ることはない）。

**Q2. 画像が多すぎて HTML が GitHub に載らない。どうするか？**
→ **縮小版と原寸版を必ず両方作り、縮小版だけ管理する。**

| | サイズ | 扱い |
|---|---|---|
| `oligo_eptA_analysis.html` | 縮小版 | **git で管理** |
| `oligo_eptA_analysis_full.html` | 原寸版 115〜120 MB | `.gitignore`。ローカルのみ |

- `6_pathway.ipynb` だけで画像 227 枚・72.7 MB あり、原寸だと **GitHub の 1ファイル
  100 MB 上限を超える**。ブラウザで開くのも実用的でない。
- 縮小版は `.qmd` の中の base64 画像を **幅 1600 px 上限・JPEG q88** に置き換えてから
  render する。**notebook 本体には手を入れない。**
- ★ **`build_html.py` は毎回2本とも作る。**片方だけ古くなる事故を防ぐため。

**Q3. 95 MB の `6_pathway.ipynb` を git でどう扱うか？**
→ **本体は `.gitignore`。写しを2つコミットする。**

| 写し | 中身 | サイズ |
|---|---|---|
| `qmd/6_pathway.qmd` | コード+markdown。GitHub で差分を読むため | 約 0.2 MB |
| `6_pathway_nofig.ipynb` | **図だけ落とした notebook。**文字出力・表・エラーは残るので数値は読める | **2.1 MB** |

- 本体をコミットすると、更新のたびに約 95 MB がリポジトリに積まれる。
  画像 227 枚（base64 で 97 MB 相当）がその大半。
- 図そのものは **HTML と `/mnt/scratch/.../pathway/`** に残るので失われない。
- `.gitignore` に入れた notebook はステージされないため、pre-commit フックが
  **明示的に拾って**両方を作り直す。図を落とす処理は `strip_figures.py`。

**Q4. GitHub を使うか／公開範囲は？**
→ **使う。Private。** 臨床データは含まないが未発表のため。アカウントは `4takeshi`。
論文発表時に Public へ切り替えられる。

**Q5. 認証は？**
→ **SSH 鍵。**`261013_Blastp_acdA` で作った `~/.ssh/id_ed25519_github` と
`~/.ssh/config` をそのまま使う（新しく作る必要はない）。

**Q6. どのファイルを管理し、どれを外すか？**
→ 「再実行で作り直せる副産物」と「GitHub の上限に触れるもの」を外した。
**scratch 側には全部そろっている**ので、外しても失われない（`CLAUDE.md` §6-10）。

| 管理する | 外す（`.gitignore`） |
|---|---|
| notebook 7本、`qmd/` の写し8本、`6_pathway_nofig.ipynb`（図なし 2.1 MB） | `6_pathway.ipynb`（95 MB） |
| `oligo_eptA_analysis.html`（縮小版） | `oligo_eptA_analysis_full.html`（原寸版） |
| `CLAUDE.md`、`README.md`、`SETUP.md`、`build_html.py`、`styles.css` | `_build/`、`.quarto/` |
| 成果物の CSV / PNG / HTML（project 直下）、`pdf/` | `annmap/`（33 MB）、`annmap_pdf/`（22 MB） |
| `refqc/`、`mappability/`、`ani/`、`logs/`、`eptA_sites/`、`ref/oligo_eptA.gtf` | `chat_preview_images/`（18 MB）、`preview_*` |
| `ai_reference/`（設計草案） | `ref/ko00001_*.keg`（4.3 MB × 3。KEGG から取り直せる） |

## 行った作業

1. **`ai_reference/`** を作り、claude.ai の設計草案2本を移した。README で現行解析との
   違い（5株→4株、パスの違い）を書いた。
2. **`build_html.py`** を作成した。
   - ★ **markdown の水平線 `---` を `***` に置き換えてから** quarto に渡す。
     quarto は markdown 中の `---` … `---` を YAML メタデータブロックとして読もうとし、
     `YAMLException: end of the stream or a document separator is expected` で落ちる（実測、42本）。
     置き換えは `_build/` に写したものだけで、**元の notebook は変えない**。
   - 縮小版・原寸版を両方作る。
3. **`styles.css`** を `261013_Blastp_acdA` からコピーした（幅の広い表が右の目次に重なる問題への対処）。
4. **`.gitignore`** を作成した（上の表のとおり）。
5. **`.githooks/pre-commit`** を作成し、`git config core.hooksPath .githooks` で有効にした。
   - notebook をコミットすると `jupytext --to qmd` で `qmd/` に写しを作り、一緒にコミットする。
   - ★ **`.gitignore` に入れた `6_pathway.ipynb` も明示的に拾う**（ステージされないため）。
   - **変換は ipynb → qmd の一方向だけ。** `jupytext --sync` は markdown セルが
     結合・消失する不具合があるため使わない（`260914_KAD_metaFL16S/SETUP.md` 参照）。
6. **`README.md`** を作成した。
7. `git init` し、最初のコミットを作成した。

## 今後の保存手順

```bash
cd /home/tanoue/261016_Xiaoxi_eptAoligo_2wFeces_Metatranscriptome

# 1. notebook を変更したら HTML を作り直す（縮小版・原寸版の両方ができる）
python3 build_html.py

# 2. 変更を記録する（notebook が含まれていれば qmd/ の写しも自動で更新される）
git add -A
git commit -m "変更内容を短く説明"

# 3. GitHub に送る
git push
```

- コミット前に `git status` で、意図しないファイルが入っていないか確認する。
- HTML を見るときは、VS Code のコマンドパレット →「Simple Browser: Show」→
  `file:///home/tanoue/261016_Xiaoxi_eptAoligo_2wFeces_Metatranscriptome/oligo_eptA_analysis.html`。
  GitHub 上では HTML として表示されない。

## GitHub への最初の接続

リポジトリ名: **`261016_Xiaoxi_eptAoligo_2wFeces_Metatranscriptome`**（Private）

### 1. GitHub のサイトで（自分で操作）

New repository → 名前に上記を入れ、**Private** を選ぶ。
**README・.gitignore・license はどれも追加しない**（空のまま作る）。

### 2. サーバで

```bash
cd /home/tanoue/261016_Xiaoxi_eptAoligo_2wFeces_Metatranscriptome
git remote add origin git@github.com:4takeshi/261016_Xiaoxi_eptAoligo_2wFeces_Metatranscriptome.git
git push -u origin main
```

SSH 鍵は `261013_Blastp_acdA` で設定済みのものがそのまま使える。
先に `ssh -T git@github.com` で `Hi 4takeshi!` と出ることを確かめる。

リポジトリ: https://github.com/4takeshi/261016_Xiaoxi_eptAoligo_2wFeces_Metatranscriptome

## 注意

- **qmd → ipynb の変換をしない。**`qmd/` は読み取り専用の写しとして扱う。
- **notebook を開いたまま外部（Claude Code など）で書き換えた場合**、VS Code が古い内容を
  表示し続けることがある。コマンドパレット →「File: Revert File」で読み込み直す。
- **`6_pathway.ipynb` は git に入っていない。**消すと復元できないので、
  scratch 側（`/mnt/scratch/.../`）にも控えを置いておくこと。
  `6_pathway_nofig.ipynb` は**図が落ちている**ので、これだけでは元に戻せない。
- **`6_pathway_nofig.ipynb` を編集しない。**元の notebook から作り直される。
