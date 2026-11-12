# eptA oligo-colonization gnotobiotic mouse metatranscriptome

4株 defined community（`Bth_VPI` / `Bth_St8` / `Bun_8492` / `Odo_St23`）＋ プラスミド `pME008` を
定着させたノトバイオートマウスの糞便メタトランスクリプトーム。
**WT 6匹 vs eptA 欠失 KO 5匹**を比べる。

- 解析の設計・決定事項・守るべき規則 → **`CLAUDE.md`**
- HTML 化と git の手順 → **`SETUP.md`**
- 全ステージをまとめた HTML → **`oligo_eptA_analysis.html`**
- リポジトリ → https://github.com/4takeshi/261016_Xiaoxi_eptAoligo_2wFeces_Metatranscriptome （Private）

## notebook（この順に実行する）

| | notebook | 中身 |
|---|---|---|
| ステージ0 | `1_genome_copy.ipynb` | ゲノムの収集 |
| ステージ1 | `2_reference_qc.ipynb` | ANI・一意配列・pangenome・mappability → **ゲート1** |
| ステージ2 | `3_build_reference.ipynb` | 統合リファレンスと GTF（eptA は3 feature に分割） |
| ステージ3-4 | `4_reads.ipynb` | fastp・bowtie2・featureCounts・QC → **ゲート2** |
| ステージ5 | `5_1_DESeq2_A_community.ipynb` | **(A) 全株プール** = 宿主が何にさらされているか |
| ステージ5 | `5_2_DESeq2_B_perstrain.ipynb` | **(B) 株ごと独立** = 各菌が自分の中で何を上げ下げしたか |
| ステージ6 | `6_pathway.ipynb` | ORA / GSEA / pathview（5ソース × 2スキーム）。**95 MB のため git では `6_pathway_nofig.ipynb`（図なし 2.1 MB）を管理** |
| ステージ6 付録 | `6_2_pathway_appendix.ipynb` | ステージ6 で図にしていなかった6経路の棒グラフとゲノム位置マップ。セル12 / セル13 と同じコードで、検定はやり直していない |
| ステージ7 | `7_genome_position.ipynb` | ゲノム座位からの解析・ゲノムブラウザ |

★ **(A) と (B) に主従は付けない。**問いが違うので両方を対等に出す。

## 主な成果物

| ファイル | 中身 |
|---|---|
| `oligo_eptA_analysis.html` | **全ステージを1つにまとめた HTML**（画像は縮小版） |
| `genome_browser_prototype.html` / `_en.html` | IGV 風ゲノムブラウザ。9,675 遺伝子・11 パネル。**locus_tag などで検索できる** |
| `panel2.html` / `panel2_en.html` | 経路レベルの2パネル対比図 |
| `eptA_samedir_DEG_loci.csv` | 同方向に3遺伝子以上連続する DEG 区間（641 行） |
| `eptA_samedir_DEG_loci_genes.csv` | その構成遺伝子（3,006 行） |
| `eptA_top10_schemeA.csv` | 株ごとの上位10区間（注釈つき） |
| `genome_annotations.csv` | **手で保守する入力。**議論で付けた領域の注釈 |
| `genome_annotations_stage456.csv` | ステージ4〜6 の議論から本データで検出したオペロン |
| `pdf/` | 図のベクタ版 |
| `ai_reference/` | claude.ai が作った設計草案（実行していない参考資料） |

## 置き場の約束

**出力はまず scratch に書き、必要なものだけ project へ複製する**（`CLAUDE.md` §6-10）。

```
/home/tanoue/261016_Xiaoxi_eptAoligo_2wFeces_Metatranscriptome/   ← PROJECT（git 側）
/mnt/scratch/tanoue/261016_Xiaoxi_eptAoligo_2wFeces_Metatranscriptome/   ← scratch（全部そろっている）
```

project にあるものは scratch にもある。例外は手で保守する `genome_annotations.csv` だけ。

## HTML の作り直し

```bash
cd /home/tanoue/261016_Xiaoxi_eptAoligo_2wFeces_Metatranscriptome
python3 build_html.py
```

★ **`git push` は自分で実行する**（`SETUP.md` 参照）。Claude Code には `git commit` までを任せる。

**縮小版と原寸版を必ず両方作る**（片方だけ古くなる事故を防ぐため）。原寸版は
`oligo_eptA_analysis_full.html` で、GitHub の 1ファイル 100 MB 上限を超えるので
`.gitignore` に入れてありローカルのみ。
