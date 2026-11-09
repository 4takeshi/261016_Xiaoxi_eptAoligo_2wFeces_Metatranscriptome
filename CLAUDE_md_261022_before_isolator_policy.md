# eptA oligo metatranscriptome — 解析設計書

プロジェクトルートに置く。Claude Code はこれを自動で読む。
**ここに書かれた判断は検討済み。再検討せず従うこと。** 変更が必要と判断した場合は、
実行せずに理由を述べて確認を取る。

対象: `261016_Xiaoxi_eptAoligo_2wFeces_Metatranscriptome`

---

## 1. 実験

gnotobiotic mouse に5株を定着、2週間後の糞便で bacterial metatranscriptome。
**eptA（lipid A phosphoethanolamine transferase、別名 lptA）が野生型の群と変異型の群を比較する。**

| 記号 | 群 | 検体 | sample name | n |
|---|---|---|---|---|
| W（= `WT`） | eptA 野生型コンソーシアム | W1–W6 | `eptAW1`–`eptAW6` | **6** |
| K（= `KO`） | eptA 変異型コンソーシアム | K1–K4, K6 | `eptAK1`–`eptAK4`, `eptAK6` | **5**（K5はサンプリングできず欠番） |

**解析での表記**（`4_reads.ipynb` セル4 で決定）:
- sample name は **`eptA` + 元の検体名**（`eptAW1` など）
- `group` 列は **`WT` / `KO`**（`W` / `K` ではない。どちらが対照か名前から自明にするため）
- DESeq2 では `relevel(group, ref="WT")` を明示する（R の factor は既定でアルファベット順に
  並ぶため、指定しないと `KO` が基準になる）

- **両群は別アイソレータで飼育。** 交雑の懸念はない
- **アイソレータ（＝飼育環境）は群と完全に交絡する。** gnotobiotic では不可避。
  → **DESeq2 の design にケージ/アイソレータ項を入れない**（群と共線）。
  論文には限界として明記する
- シーケンス: **150 bp PE**、MGI/DNBSEQ、Azenta。
  **アダプター未除去**（P7 read1 `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA` /
  P5 read2 `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT`）
- 深度は 43.6M〜90.8M read pairs（2.1倍の開き）

## 2. 構成株と変異の構造（配列から確認済み）

| 株ID | 株 | 作製 | 変異法 | 構造 |
|---|---|---|---|---|
| `Bth_VPI` | *B. thetaiotaomicron* VPI-5482 | Xiaoxi | 二重交叉 | **触媒ドメインのみ 503 bp 欠失**。カセットなし |
| `Bth_St8` | *B. thetaiotaomicron* CE91-St8 | 田之上 | 二重交叉 | **ほぼ全長 1,690 bp 欠失**（5'端数十塩基は前遺伝子と重複のため残存）。カセットなし |
| `Bun_8492` | *B. uniformis* ATCC 8492 | Xiaoxi | 二重交叉 | **触媒ドメインのみ 575 bp 欠失**。カセットなし |
| `Odo_St23` | *Odoribacteraceae* CE91-St23 | 田之上 | CAST | **1,445 bp 挿入**（150bp Left Flank + tetX1 + PcfxA + 57bp Right End） |
| `Odo_St24` | *Odoribacteraceae* CE91-St24 | 田之上 | CAST | 同一の挿入。**ゲノムは St23 とほぼ同一（下記）** |

**プラスミド pME008（= pSL5053、15,808 bp、Type I-F CAST）は *Odoribacter* のみが持つ。**
- W群の Odo: エンプティーベクター（スペーサーなし）
- K群の Odo: eptA標的スペーサーを持つターゲティングベクター
- **両方とも tetX1 と CAST機構一式（Cas6/7/8, TniQ, TnsA/B/C）を保持**
- → **tetX1 は群の判別に使えない**
- → ~~**ただしスペーサー領域の 32 bp は判別に使える**（2026-10-20 確定。下記）~~
  **★ RNA-seq では使えない（2026-10-22 確定。下記）。** 配列としては判別できるが、
  **転写されないのでリードが来ない**（深度 1.7、p = 0.85）。群の判別は eptA `_deleted` で行う

### pME008 のスペーサー領域（~~★群の判別子~~ → **RNA-seq では使えない。2026-10-22 確定**）

> **★ 訂正（2026-10-22、`4_reads.ipynb` セル14-7）**
> **この領域は RNA-seq の判別子にならない。** 全11検体で **深度 1.7**、W/K = **0.8倍**、
> **p = 0.85**。転写されない領域は RNA-seq では見えないので、カバレッジ差が出ない。
> **→ 群の判別は eptA の `_deleted` feature で行う**（§4。W/K = 95〜4,191倍、
> Salmon TPM では K群 全5検体が 0.00）。
> 以下の記述は**配列設計としては正しい**ので、DNA レベルの確認には使える。

W群と K群の**唯一の配列差**。`pME008.gb` 上の座標と配列を確定済み
（`ref/pME008_spacer.tsv`、`3_build_reference.ipynb` セル4-2）。

```
pME008 (pSL5053, 15,808 bp) 上の 9,421–9,452（+鎖、32 bp）
  W群 empty  : ATAAC AGTCTTCCTGTCTTGTCAGCTAGGTGAAGACT GTGAA
  K群 mutant : ATAAC GGTACTGATCTTACACGCCTTCCAGCTCGTTT GTGAA
              └5bp┘ └────────── 32 bp ───────────┘└5bp┘
```

- **挿入ではなく置換**（32 bp ↔ 32 bp）。プラスミド全長の差 1 bp は環状配列の線状化による
- GenBank に `replace_spacer` という feature が実在し、前後を `V. cholarae repeat`
  （CRISPRリピート）に挟まれた repeat–spacer–repeat 構造。内部の BbsI 2か所が
  golden-gate 部位で、W群の 32 bp は stuffer（`GAAGAC` を含む）
- リファレンスには**エンプティー側を入れる**ので、**K群ではこの領域に重なるリードが
  マップされずカバレッジが落ちる**（150 bp リードなので前後を含む約180 bp の窓）
- 染色体の CAST 挿入ジャンクション検出より単純（カバレッジを見るだけ）。
  *Odoribacter* の遺伝子型を直接確認できる
- **Bacteroides 3株（VPI・St8・B.uni）はプラスミドを持たない。** 二重交叉でマーカーレス。
  St8 の WT 配列ファイルにある `spacer_st08_eptA1/2` は設計時の注釈であり、実体はない（確認済み）
- → **pME008 にマップされたリードは Odoribacter のプラスミド発現**として解釈する

## 3. ゲート1の結果（実施済み・確定）

| ペア | ANI | 1kb以上の固有配列（少ない方） | 判定 |
|---|---|---|---|
| Bth_VPI ↔ Bth_St8 | 98.03% | **199,018 bp**（206 CDS） | GO |
| Odo_St23 ↔ Odo_St24 | 99.99% | **0 bp** | **STOP** |

**Odo_St23 と St24 は同一株とみなす。**
- リファレンスには **St23 のみ**を入れる。St24 は入れない
- カウントは「*Odoribacter* sp.（CE91-St23/St24）」という1単位。
  論文・図でもそう表記し、2株を区別していないことを明示する

**Bth ペアは株レベルで解析可能。mappability の実測で確定済み（2026-10-20）。**

| 株 | 全CDS | 区別可能（mappability ≥ 0.5） | 固有領域のみなら |
|---|---|---|---|
| Bth_VPI | 4,803 | **2,369 (49.3%)** | 206 |
| Bth_St8 | 4,797 | **2,337 (48.7%)** | 230 |
| Bun_8492 | 3,690 | 3,613 (97.9%) | — |
| Odo_St23 | 4,620 | 4,389 (95.0%) | — |

**「固有領域の206 CDSに限定」ではなく、その10倍以上が使える。**
共通領域でも、2% の塩基分岐が 150 bp 中に十分な数入っている CDS は区別できる。

分布は**完全な二峰性**（Bth 2株は 0 付近と 1.0 付近に山、中間は少数）なので、
閾値を 0.5 に置いても 0.9 に置いてもほぼ同じ集合が選ばれる。閾値が結論を左右しない。

測定は本番リファレンス（pME008 + マウス decoy 込み、2.76 GB）に対して実施。
**マウス decoy は細菌CDSの判別に影響していない**（ステージ1の細菌4株のみの測定との
平均差 0.0002、`mouse` と混同される細菌CDSは0件）。
**pME008 の 13 CDS はすべて mappability 1.000** なので、pME008 にマップされたリードは
確実にプラスミド由来と解釈できる。

値は GTF の `mappability` 属性に入っている（`3_build_reference.ipynb` セル23・26-28）。

**タンパク質レベルで検証済み（2026-10-19）**。固有領域のCDSを St8 プロテオームに
diamond blastp したところ、ゲノム全体の 3,527/4,645 CDS (75.9%) にオーソログがあるのに対し、
固有領域は R1 0/64・R2 0/51・R3 0/34・R5 0/1・p5482 0/38（PUL のみ 1/7）。
**76% 対 0% の対比により、これらは真の固有領域と確定。199,018 bp を減額しない。**

> 当初「固有領域に 30S ribosomal protein S17 があるのはおかしい」としたのは誤り。
> VPI は S17 を2コピー持ち、必須コピー（WP_008762040.1 / BT_RS13740、3.37 Mb の
> リボソームタンパク質オペロン内）は St8 の GKH19591.1 と100%一致する。固有領域にあったのは
> 同一性 50.8% のパラログ（WP_055235409.1 / BT_RS00025）。「必須単一コピー」という前提が
> この株では成り立っていなかった。

**R2（68,430–111,285、42.9 kb、51 CDS が全て固有）は生物学的に注目に値する。**
RagB/SusD + ABC transporter + DDE transposase を含み、contig端とは無関係の位置にある。
R4 の α-L-fucosidase を含む PUL と合わせて、両株の多糖利用能の差の候補。

## 4. ★ eptA は DE テーブルから読めない

**4つの変異がすべて異なる挙動をし、どれもカウント0にはならない。**

**座標は GFF の CDS 基準。** 断片（SnapGene/GenBank）由来の座標ではない。
マッピングは GFF の CDS 全長に対して行われるので、カウントも CDS 基準になる。

| 株 | eptA の locus_tag | 鎖 | CDS長 | K群での eptA カウントの挙動 |
|---|---|---|---|---|
| Bth_VPI | `BT_RS09185` | **−** | 1,704 bp | 504 bp 欠失 → 残り **1,200 bp**（70.4%、log2FC ≈ **−0.51**） |
| Bth_St8 | `CE91St8_20540` | + | 1,704 bp | 1,689 bp 欠失 → 残り **15 bp** のみ（log2FC ≈ **−6.83**）。実質ゼロだが完全な0にはならない |
| Bun_8492 | `NQ510_RS14730` | + | 1,701 bp | 575 bp 欠失 → 残り **1,126 bp**（66.2%、log2FC ≈ **−0.60**） |
| Odo_St23 | `CE91St23_40970` | **−** | 1,761 bp | 挿入点は **CDS 337 位の直後**。5'側 337 bp は残り、3'側 1,424 bp は転写終結次第（log2FC 0 〜 **−2.39**） |

> **Bth_VPI の数値が当初の 1,002 bp / −0.59 から変わっている**（2026-10-20）。
> SnapGene の `lptA` feature（断片上 1..1506）が部分的な注釈で、
> GFF の実際の CDS は 1,704 bp だった（198 bp の差。マイナス鎖なので遺伝子の3'末端側）。
> `3_build_reference.ipynb` セル21 で食い違いを検出・表示する。
>
> **鎖は株ごとに違う**（VPI と Odo_St23 がマイナス鎖）。挿入の 5'/3' 分割で
> これを間違えると、どちらが転写終結の影響を受ける側かが逆になる。
>
> **DDBJ の2株（Bth_St8 / Odo_St23）は product が `membrane protein`** で
> `phosphoethanolamine` を含まないため、product 名では eptA を特定できない。座標で特定する。

**そのまま DE 結果に載せると「中程度に低下した遺伝子」と誤読される。** 遺伝子型による
見かけの低下であって、発現制御の変化ではない。

### 対処: GTF に3つの feature を置く

`gene_id` は `<元の locus_tag>_<接尾辞>` にする（`3_build_reference.ipynb` セル21）。

```
<locus_tag>            元のCDS全長            → DE解析から除外する（exclude_from_DE=true）
<locus_tag>_retained   両群に共通する領域      → ★ここで比較する。差が出れば本物の制御変化
<locus_tag>_deleted    K群で欠失している領域   → 遺伝子型の確認用（K群でゼロになるはず）
```

**★ 残存領域が2区間に分かれても、どちらも同じ `gene_id`（`_retained`）にする。**
欠失が CDS の内部にあると残存は 5'側 と 3'側 の2区間になるが（Bth_VPI・Bun_8492 が該当）、
featureCounts は同じ gene_id の行を合算するので、1つのカウントとして扱える。

Odoribacter は挿入なので2分割する。

```
<locus_tag>_5p    挿入点より5'側   337 bp  → 両群で転写される
<locus_tag>_3p    挿入点より3'側 1,424 bp  → K群でゼロなら挿入による転写終結、
                                             残っていれば読み抜け。これ自体が情報
```

3'側が CDS の 81% を占めるので、**ここが K群で落ちるかどうかが重要な観察項目。**

座標は確定済み（2026-10-20）。`refqc/eptA_mutation_spec.tsv` と
`3_build_reference.ipynb` セル7・9・10・21 を参照。

## 5. リファレンス構成

```
ref/oligo_eptA.fna
  ├── Bth_VPI~<contig>      野生型ゲノム（NCBI）
  ├── Bth_St8~<contig>      野生型ゲノム（NCBI）
  ├── Bun_8492~<contig>     野生型ゲノム（NCBI）
  ├── Odo_St23~<contig>     野生型ゲノム（NCBI）   ※St24 は入れない
  ├── pME008~pSL5053        15,808 bp            ※Odoribacter が両群で保持。入れないと unmapped
  └── mouse~<chr>           decoy
```

- **両群とも同じリファレンスを使う。** 別々に作ると遺伝子IDが対応せず比較できない
- **配列を連結して1本にしない。** ファイルを連結した multi-FASTA にし、contig は独立レコードのまま
- **contig名の先頭に株IDを付ける**（`株ID~contig名`）。株間で contig 名が衝突するのを防ぐ
  （Bakta を使う場合は `contig_1` 等に改名されるため特に必要）
- **マウスゲノムを decoy として追加する。** 宿主RNAが細菌遺伝子へ誤マップするのを防ぐ
- pME008 は**スペーサー 32 bp が置換されているだけ**なので、エンプティーベクター側を
  1本入れればよい（§2）。GTF には `CAST_spacer`（9,421–9,452）を feature として入れ、
  ステージ4でこの領域のカバレッジを群ごとに比較する

## 6. 絶対に守ること

1. **`featureCounts` は必ず `-s 2`（reverse-stranded）。** Rsubread なら `strandSpecific=2`。
   同系統の過去データで rev/fwd = 20〜55:1 を実測済み。未指定だとアンチセンスが約4%誤カウントされる
2. **DESeq2 に入れるのは生カウントのみ。** 係数で割って `round()` してから渡さない。
   補正が必要なら `normalizationFactors()` に offset として渡す。
   過去プロジェクトではこの誤りで検定可能遺伝子の30%と DEG の40%を失っていた
3. **CDS のみで DE 解析する**（GTF の `biotype == "CDS"`）。rRNA/tRNA/ncRNA は**消さずに
   GTF に残して数え、DE からだけ外す**。過去データでは rRNA が assigned read の 37〜56% を占めた。
   rRNA feature が無いと ①除去効率を報告できない ②rRNA リードが隣接CDSに漏れる
   ③size factor が歪む
4. **`enricher()` には必ず `universe` を渡す。** 背景は「検定された遺伝子（padj が NA でない）
   のうち対象集合に属するもの」。渡さないと「変動した機能」ではなく「その株が持つ機能」が出る
5. **重複除去（Picard MarkDuplicates 等）をしない。** RNA-seq では高発現遺伝子が
   正当に重複リードを生むため、入れると高発現遺伝子だけが壊れる
6. **SortMeRNA を使わない。** 全部マップしてから CDS に絞る。rRNA% が QC 指標として自動的に得られる
7. **eptA を DE 結果にそのまま載せない**（§4）

## 7. ステージ

### ステージ2 — リファレンス構築（**完了 2026-10-20**。`3_build_reference.ipynb`）

eptA座標確定（minimap2）→ CheckM2 → **既存アノテーション（NCBI/DDBJ の GFF）**
→ eggNOG-mapper → contig改名 → 統合FASTA + mouse decoy
→ 統合GTF（株ID・biotype・KO・**eptAの3分割feature**・mappability を属性に持たせる）

**できあがったもの**（`/mnt/scratch/.../261016_Xiaoxi_eptAoligo_2wFeces_Metatranscriptome/`）

| パス | 内容 |
|---|---|
| `ref/oligo_eptA.fna` | 2.76 GB / 74 contig（細菌4株 12 + pME008 1 + マウス 61） |
| `ref/oligo_eptA.gtf` | 10.3 MB / 18,258 行。**DE 対象 17,922**（CDS かつ `exclude_from_DE=false`） |
| `ref/oligo_eptA_cds.fna` | 転写産物 17,919 本 |
| `ref/bt2/` | bowtie2 index 3.80 GB（contig 74 で検算済み） |
| `ref/salmon/` | decoy-aware index 6.7 GB。**定量される転写産物は 17,182 本**（= CDS 17,919 − 重複 737。`quant.sf` の行数で検算済み 2026-10-22）。`info.json` の `num_refs` 17,256 は index 内部の総数で転写産物数ではない |
| `ref/pME008_spacer.tsv` | スペーサー領域の座標（9,421–9,452）。**RNA-seq では判別子にならない**（§2 の訂正。DNA レベルの確認用） |
| `refqc/checkm2/` | 4株とも Completeness 99.95%以上 / Contamination <5% |
| `refqc/eptA_mutation_spec.tsv` | eptA 変異のゲノム座標（検算済み） |
| `eggnog/all4.emapper.annotations` | 17,341 件。COG 85〜92% / KO 37〜44% |
| `mappability_full/cds_mappability.tsv` | 本番リファレンスでの mappability |

**featureCounts では `-t transcript -g gene_id -s 2` を指定する**（feature 型が `transcript`）。

**★ Salmon を使うときの注意。** 配列が完全に同一の転写産物 **737 本**が除外されている
（`ref/salmon/duplicate_clusters.tsv`）。内訳は **Bth_St8 ↔ Bth_VPI 463 件**、
Odo_St23 自株内 88 件、**Bth_St8 ↔ Odo_St23 67 件**など。
除外された側には**カウントが付かない**ので、株ごとの比較で「片方の株だけ発現0」が出たら
まずこのファイルを確認する。

> **★ 実測（2026-10-22、セル16-5）: 除外は `Bth_VPI` に偏っている。**
> 737本の株別内訳は **`Bth_VPI` 504 / `Odo_St23` 155 / `Bth_St8` 47 / `Bun_8492` 31**。
> `minMQS=0` のカウントで見ると**除外分は `Bth_VPI` の 10.3%、`Bth_St8` の 0.2%**で、
> **これだけで `VPI:St8` が 0.659 → 0.593 に下がる。**
> **→ Salmon の `VPI:St8` 絶対値を根拠に議論してはいけない。群間比較だけが意味を持つ。**
> （なお `--fraction -M` で 1/2 ずつ配分されるという想定は誤り。bowtie2 が `NH` タグを
> 出さないため featureCounts の按分は動かない。§7 ステージ3 の訂正を参照）

**両者の相関を取って乖離を見る**（下記ステージ3）のが重要。

**アノテーションは既存の locus_tag を使う。Bakta は欠けている株だけ。**
`BT_RS13740` のような既存IDは文献で引けるが、Bakta の `Bth_VPI_00135` はどこにも無い。
*B. theta* の PUL 文献との照合でこの差が効く。NCBI/DDBJ の GFF には rRNA・tRNA も入っている。
eggNOG の入力は **`_translated_cds.faa`**（`_protein.faa` は同一配列を1レコードに統合するため
パラログを取りこぼす）。

### ステージ3 — サンプル処理（**完了 2026-10-22**。`4_reads.ipynb`）
FastQC/MultiQC → fastp（アダプター明示、`--length_required 50`）
→ Bowtie2（`-X 1000`）→ samtools → featureCounts（`-t transcript -g gene_id -s 2 -p --countReadPairs -B`、`minMQS` 2通り）
→ **Salmon を並走**（selective alignment、decoy に全ゲノム、`-l ISR`）。Bth 2株が近縁なため必須

> **★ 訂正（2026-10-22、セル15）: featureCounts の `--fraction -M` は効かない。**
> featureCounts は多重マップを **`NH` タグ**で判定するが、**bowtie2 は既定出力に `NH` を付けない**
> （`AS` / `XS` / `NM` / `MD` のみ）。全リードがユニーク扱いになり、`-M` も `--fraction` も何も起こさない。
> **→ 多重マップの扱いは次の3通りで揃える。**
>
> | 見方 | 実装 | 用途 |
> |---|---|---|
> | 捨てる | featureCounts `minMQS = 10` | **ステージ5 の主解析** |
> | 最良ヒットに全部入れる | featureCounts `minMQS = 0` | 比較用 |
> | **確率的に按分** | **Salmon の EM** | **感度解析** |
>
> **実測（セル16-5）**: `Salmon / featureCounts(minMQS=10)` の比の中央値は
> `Bth_St8` **2.161** / `Bth_VPI` **1.708** に対し `Bun_8492` 0.835 / `Odo_St23` 0.855。
> **EM が拾っているのは近縁2株の共通領域だけ。**

**実測の要点（11検体）**: fastp 残存率 99.35〜99.82% / bowtie2 アライン率 98.25〜98.96% /
rRNA% 0.10〜0.29（`minMQS=10`）・2.40〜12.65（`minMQS=0`）/ Salmon CDS マップ率 34.5〜46.1% /
size factor 0.65〜2.36（3.6倍。目安の5倍以内）。**所要時間**: bowtie2 223分 / featureCounts 36分 / **Salmon 11分**。

### ステージ4 — QCと遺伝子型確認（**ゲート2 通過。2026-10-22**。`4_reads.ipynb` セル17）

**10項目中 通過8 / 使用不可1 / 未実施1。**

| 項目 | 期待 | **実測（2026-10-22）** | 判定 |
|---|---|---|---|
| **eptA 欠失区間のカバレッジ**（VPI/St8/B.uni） | W群 = あり / **K群 = ゼロ** | 深度比 95〜4,191倍（p=0.0043〜0.0078）。**Salmon TPM は K群 全5検体 0.00** | **通過** |
| **eptA 残存区間のカバレッジ** | 両群にあり。ここで比較する | 両群にあり。欠失境界を補正した比は VPI 1.96〜2.38 / Bun 0.66〜1.90 | **通過** |
| **Odo 染色体のCAST挿入ジャンクション** | W群 = なし / **K群 = あり** | `_3p` が W/K = 218倍（p=0.0043）。`_5p` は両群にあり（挿入で転写が止まっている） | **通過** |
| **pME008 スペーサー領域（9,421–9,452）のカバレッジ** | W群 = あり / K群 = ゼロ | **深度 1.7、W/K = 0.8倍、p = 0.85** | **★ 使用不可**（§2 の訂正。eptA `_deleted` で代替） |
| **トランスポゾン–染色体ジャンクションの全ゲノム検索** | W群 = なし | **挿入点は転写されないので RNA-seq では検出できない。DNA（WGS）が必要** | **★ 未実施**（ステージ3 の守備範囲外） |
| tetX1 | **両群にあり。判別に使えない** | 全11検体で検出（768〜3,572 リード / TPM 17.1〜46.0） | 通過（判別子としては使わない） |
| pME008 由来リードの割合（群ごと） | 群間でずれていないか | 全検体 0.006〜0.016%。size factor に影響する量ではない | **通過** |
| 株別CDSリード数・シェア・unique率 | 下記の除外基準 | 4株とも基準の3〜60倍の余裕 | **通過** |
| **Bth_VPI : Bth_St8 の存在比（群ごと）** | 大きく変わっていないか | **4つの見方すべてで群間差なし（p = 0.54〜0.79）。** 検出力は相対差40〜50%まで | **通過**（検出力の制約つき） |
| サンプル別 mapping率 / rRNA% / multimap率 | — | アライン率 98.25〜98.96% / rRNA% 0.10〜0.29（`minMQS=10`）/ **decoy行き率に群間差なし（p=0.66）** | **通過** |

**★ 群集組成には群間差がある**（Bth WT 42.0% vs KO 32.0%、Bun WT 54.8% vs KO 65.1%、ともに p = 0.0043）。
**アイソレータと交絡している（§1）ので eptA の効果とは言えない。** ステージ5・6 の解釈に持ち越す。

**除外基準: 総CDSリード 200,000 未満、または全CDSリードの 0.1% 未満の株は株レベル解析から外す。**
過去プロジェクトでは全8検体で1リードしかない株が最後まで解析に残り、論文の主張を支えていた。

### ステージ5 — 統計（DESeq2 を2本）

```
(A) 全株プール DESeq2   → 群集としての機能出力。転写産物量の変化を含む  ＝「総体」
(B) 株ごと独立 DESeq2   → 株内の発現制御。転写産物量の変化を含まない    ＝「細菌ごと」
```

**★ (A) と (B) に主従はつけない（2026-10-22 決定）。見ている問いが違うので、両方を対等に出す。**

| | **(A) 全株プール** | **(B) 株ごと独立** |
|---|---|---|
| 問い | **群集として何の転写が増えたか** | **各菌が自分の中で何を上げ下げしたか** |
| size factor | ライブラリ全体で1組 | その株の中だけで各1組 |
| 株のシェア変化 | **含む** | **含まない** |

**最も情報量が多いのは「片方のパネルでのみ有意」な経路。**
過去プロジェクト T19 の例: **PTS は (A) padj 1.1e-12 / (B) padj 0.10** —
実体は *Enterococcus hirae* のシェアが 0.56%→1.75% と3倍になっただけ。
逆に **RNAポリメラーゼは (B) padj 6.3e-05 / (A) n.s.** — 菌量と無関係に各菌が増殖装置を絞っている。

**※ 用語**: ここでいう「シェア」はメタトランスクリプトームの**転写産物量**であって菌体数ではない。

**★ ステージ3 からの申し送り（2026-10-22）**

| # | 内容 |
|---|---|
| 1 | **入力は `counts/featureCounts_q10.csv` が主解析。** 感度解析は `salmon_quant/*/quant.sf` を **tximport** で読む（`EffectiveLength` の長さオフセットが要る。`counts/salmon_*.csv` は比較用） |
| 2 | **DE 対象から外す**: 遺伝子型マーカー（`_deleted` / `_3p` / `_5p` / `_retained`）。**`BT_RS09185_retained` も**（KO で 2.2倍低いが交絡と区別できない） |
| 3 | **size factor は CDS のみで計算する**（rRNA・tRNA・ncRNA を含めない） |
| 4 | **`eptAW4` を除いた感度解析**（§8） |
| 5 | **DE 結果に「CPS 座位」の列を付ける**（§8。セル15-10-6 で広げた範囲を使う） |
| 6 | **`product` だけでなく `Preferred_name` / `PFAMs` / `Description` を見る**（`product` だけでは `neuA` も `fnlA` も見逃す） |
| 7 | **Salmon の `VPI:St8` 絶対値は使わない。** 重複クラスタ737本が `Bth_VPI` に偏り（`Bth_VPI` の 10.3% / `Bth_St8` の 0.2%）、それだけで比が10%下がる |
| 8 | **Salmon の `--meta` 版は不要**（群間比較では VBEM と完全に同じ p。保存はしてある） |

- **★ 技術的な注意: (A) では低存在量株の p 値が出やすい。**
  DESeq2 は dispersion を全遺伝子の trend から借りるので、(A) では**低存在量株
  （今回は `Odo_St23`、3.8%）の dispersion が高存在量株（`Bun_8492`、72.4%）に引っ張られて
  過小評価されうる**（＝偽陽性が出やすい）。
  **これは (A) を格下げする理由ではなく、(A) の `Odo_St23` の結果に注記を付ける理由。**
  (B) では株の中で dispersion を推定するのでこの問題は起きない
- **(A) に株ごとの offset を足す方式は採らない。** 株全体が丸ごとシフトして全経路が汚染される
- `design = ~ group`（アイソレータ項は入れない）
- `lfcShrink(type="apeglm")`。FCでランク付けする箇所は必ず収縮後の値を使う
- prefiltering: 最小群サイズ以上のサンプルで count >= 10
- n = 6 vs 5 の不均衡は DESeq2 が扱える

### ステージ6 — 機能解析

#### ★ 最終成果物の形（2026-10-22 決定）

過去プロジェクト T19 の成果物 **`LP食での経路シフト`**（2026-10-16）と同形式の **2パネル図**。

```
            経路名        │  (A) 群集としての機能出力  │  (B) 株内の発現制御
  ────────────────────────┼───────────────────────────┼──────────────────────
  鞭毛形成  ko02040 n=187 │        ●───→ +0.161       │      ●─────→ +0.306
  PTS糖取込 ko02060 n=586 │        ●───→ +0.186 ★     │      ○→ +0.053  n.s.
```

- **左右で x 軸・行順を共通**にして、同じ経路が横に並ぶようにする
- **スコア = 2×(AUC−0.5)**。経路内の遺伝子の **Wald 統計量**が、それ以外の遺伝子と比べて
  どちらに寄っているか。+1 で完全に上側、−1 で完全に下側
- **片方のパネルでのみ有意な経路をハイライトする**（最も情報量が多い）
- ホバーで**担い手の株**を出す。**順位は遺伝子数ではなく転写量（正味リード増減）**

**→ ステージ5 で保存しておくもの（これが無いと図が作れない）**

| 保存するもの | 用途 |
|---|---|
| **遺伝子ごとの Wald 統計量**（(A)(B) 両方） | パネルのスコア計算 |
| 遺伝子ごとの LFC・padj（(A)(B) 両方） | 表・担い手の判定 |
| **群ごとの正規化カウント平均**（(A)(B) 両方） | **担い手の順位＝正味リード増減** |
| 遺伝子 → 株 / 経路 / CPS座位 の対応 | 行の分解 |

**※ (B) の p 値は、株ごとに DESeq2 を回したあと p 値をプールして BH 補正を1回かける**
（T19 の成果物と同じ作法）。

**★ 今回と T19 の決定的な違い**: T19 は食餌介入（AIN93G vs LP）で交絡が無かったが、
**今回は群とアイソレータが交絡している（§1）。左パネル（群集としての機能出力）は
その交絡をそのまま含む。** 図の注記に必ず書くこと。

#### 手法

- **fgsea が主**（Wald統計量でランク付け）。ORA は従で `universe` 必須
- KEGG pathway に加えて **KEGG module・CAZy/dbCAN（PUL）・GMM/GBM・TCDB**。
  *Bacteroides* 主体の群集なので多糖利用が主役になりうる
- global map（ko01100/01110/01120/01200/01210/01212/01230 等）と
  eggNOG が細菌に誤付与する真核・疾患 pathway を除外
- **KEGG のラベルを鵜呑みにしない**。ko00195「Photosynthesis」の実体は F0F1 ATP synthase、
  ko00120/00121「bile acid biosynthesis」の実体は BSH。図には実体を書く
- 有意経路はゲノム座標順に LFC をプロットし、1オペロンが効いているだけでないか確認
- **担い手の株は遺伝子数ではなく転写量（正味リード増減）で順位付ける。**
  少数の高発現遺伝子を持つ株は遺伝子数では埋もれる
- **★ CPS 座位をそのまま経路解析に入れない**（§8。2026-10-22）
  `Bth_St8` の CPS 座位の切り替わりは**実体が「1回の相転換」なのに遺伝子数としては数十個**
  になるので、「莢膜多糖生合成」「糖ヌクレオチド代謝」が **fgsea でも ORA でも必ず上位に来て
  過大評価される。** **座位をまとめて1事象として扱うか、除いた場合の結果も併せて出す。**
  図・本文では「**CPS 座位の切り替わり（アイソレータと交絡）**」と明示する

## 8. 未決事項

- [x] ~~eptA の CDS 座標と欠失/挿入区間の照合~~ → 完了。`eptA_mutation_spec.tsv`。
      残り: Bth_VPI と Bun_8492 のゲノム座標（ステージ2 §C の minimap2 で確定）
- [x] ~~mappability の実測~~ → 2026-10-20 完了。本番リファレンス（pME008 + マウス decoy 込み）で
      測定し、GTF の `mappability` 属性に反映済み。**Bth 2株で約2,350 CDS が使える**（§3）
- [x] ~~Bth_VPI 固有領域のタンパク質レベル確認~~ → 2026-10-19 完了。真の固有領域と確定（§3）
- [x] ~~Bakta DB の方針決定~~ → 2026-10-20。**既存アノテーションを使う。Bakta は欠けた株だけ**
- [x] ~~Bth_St8 も pME008 を持つか~~ → 2026-10-20。**持たない。** pME008 は Odoribacter のみ
- [x] ~~`eptA mutant sites` フォルダのサーバ転送~~ → 2026-10-20 完了。`eptA_sites/`
- [x] ~~各株のアセンブリ accession の確認~~ → 2026-10-20 完了（§11）
- [x] ~~rRNA 除去効率が群間でずれていないか~~ → 2026-10-22 否定。decoy行き率 WT 54.90% vs
      KO 56.91%、**p = 0.66**（`4_reads.ipynb` セル16-6）
- [ ] **★ CPS 座位の切り替わり（最重要の未決事項。2026-10-22 発見）**
      PCA の PC1 の正体は、`Bth_St8` の**莢膜多糖（CPS）生合成座位が群間で入れ替わっている**こと
      （`4_reads.ipynb` セル15-10-3〜15-10-8）。
      3座位すべてが `NusG` + `UpxZ` ドメイン対（= `upxY`/`upxZ`）を持ち、排他的に使われている。
      - `CE91St8_26020`–`26140`（13.2 kb、L-ラムノース + O-アセチル化）**W群で高い**
      - `CE91St8_29520`–`29880`（35.4 kb、シアル酸型ノヌロソン酸 + アミノ糖 + グルクロン酸）**K群で高い**
      - `CE91St8_37520`–`37680`（18.9 kb、L-FucNAc + グルクロン酸）**W群で高い**

      **eptA の効果とは言えない**（①*Bacteroides* の CPS は確率的に相転換する、
      ②アイソレータと交絡、③「1遺伝子ずつ」でなく「座位ごと」の変化）。
      **検証には同一アイソレータ内で WT/KO を飼育するか、単培養での比較が要る。**
      **ステージ5・6 では DE と経路解析を強く歪めるので、座位の列を付けて別扱いにする。**
- [ ] **`eptAW4` を外れ値として扱うか**（2026-10-22）
      **3つの独立な指標**で浮かぶ: RLE IQR 1.37（他は 0.49〜1.01、セル15-8-2）、
      PCA の PC2/PC3（セル15-9）、Salmon の整合率 95.66%・ISR率 95.69%（ともに最低値、セル16-2）。
      **除外はしない。ステージ5 で除いた感度解析を必ず行う。**
- [ ] **`VPI:St8` の検出力不足**（2026-10-22）
      群間差なしだが、n=6 vs 5 で**検出できるのは相対差 40〜50% 以上**。小さい差は見えていない
- [ ] **アノテーションが株で非対称**（2026-10-22、セル15-4）
      RefSeq 株（VPI / Bun）には `RNase_P_RNA` / `riboswitch` / `pseudogene` があるが、
      GenBank 株（St8 / St23）には無い。**株間で biotype 別の比較をしない**

## 9. 報告の作法

各ステージの終わりに、出した数字（表で）、判断が必要な箇所、想定と違った点を簡潔に報告する。
**ゲート（ステージ1終了時、ステージ4終了時）では必ず止まり、確認を取ってから次に進む。**

### notebook の作法（`260914_KAD_metaFL16S/CLAUDE.md` と共通）

- **一度の指示につき、必ず1セルのみを追加・編集する。** 複数セルをまとめて作らない
- 明示的に「次に進んで」と言われるまで、次の解析ステップを提案・追加しない
- **セルは自分で実行しない。** 実行と結果確認はユーザーが行う
- 新しいセルを追加する前に、直前のセルの出力を確認し、意図通りかどうか一言コメントする
- **「これは残す」と決まった結果セルの直後には、考察を書く Markdown セルを追加する**
- 新しく作るセルは、そのセル単体でも実行できるよう、依存する前段を省略せず含める
- コードには `#` で「何をしているか」「なぜそうするか」を書く。試行錯誤が多い段階なので、
  後から読み取れる形を優先する
- **表は `show(df)`（HTML テーブル）で出す。** `print(df.to_string())` は日本語列名で桁が崩れる

---

## 10. 実行環境（hondaserver3）

旧解析（33株、2023–2024年）は **hondaserver2** で実施。本プロジェクトは hondaserver3 で
新規に構築した環境を使う。構築記録は `/home/tanoue/261016_metatranscriptome_analysis_enviroment/`。

| | 内容 |
|---|---|
| conda 環境 `metatx` | fastp 1.4.0 / FastQC 0.13.0 / MultiQC 1.35 / SortMeRNA 7.1.0 / STAR 2.7.11b / featureCounts 2.1.1 / Salmon 2.8.0 / samtools 1.24 / bowtie2 2.5.5 / skani 0.4.0 / SeqKit 2.14.0 / fastANI 1.34 / minimap2 2.31 |
| conda 環境 `checkm2` | CheckM2 1.1.0（DB: `/db/checkm2/CheckM2_database/uniref100.KO.1.dmnd`） |
| conda 環境 `eggnog-mapper` | emapper 2.1.12 / eggNOG DB 5.0.2（`EGGNOG_DATA_DIR=/db/eggnog-mapper`） |
| R 4.4.1 / Bioc 3.20 | DESeq2 1.46.0 / apeglm / ashr / Rsubread / tximport / edgeR / limma / tidyverse |
| **R 専用ライブラリ** | `/home/tanoue/R/metatx-lib-4.4`：ggplot2 3.5.2 / igraph 2.3.4 / fgsea 1.32.4 / clusterProfiler 4.14.6 / enrichplot 1.26.6 |

- **R notebook の冒頭で必ず実行する**:
  `.libPaths(c("/home/tanoue/R/metatx-lib-4.4", .libPaths()))`
  user lib の ggplot2 は **gggenomes 用に 3.5.1 で固定**しているため分離してある。
  gggenomes を使う notebook ではこの行を書かない（同一 R セッション内での混在は不可）
- notebook から `conda activate` は効かない。環境の `bin` を絶対パスで呼ぶか `mamba run -n <env>`
- Bakta 1.8.1 が base 環境にあるが **DB 未設定**。§7 のとおり既存アノテーションを使うので触らない

## 11. ファイルの所在

```
/home/tanoue/261016_Xiaoxi_eptAoligo_2wFeces_Metatranscriptome/   ← PROJECT（notebook・小さい表。git 側）
  ├── CLAUDE.md
  ├── 1_genome_copy.ipynb        ステージ0: ゲノム収集
  ├── 2_reference_qc.ipynb       ステージ1: ANI・一意配列・ゲート1（mappability は未実施）
  ├── 3_build_reference.ipynb    ステージ2
  ├── 4_reads.ipynb              ステージ3（完了 2026-10-22。206セル）
  ├── 5_1_DESeq2_A_community.ipynb   ステージ5 (A) 全株プール = 群集としての機能出力（総体）
  ├── 5_2_DESeq2_B_perstrain.ipynb   ステージ5 (B) 株ごと独立 = 株内の発現制御（別々）
  ├── 6_pathway.ipynb            ステージ6（(A)(B) 2パネルの経路図）
  ├── ani/                       ステージ1の結果（scratch から複製）
  ├── eptA_sites/                eptA WT/変異体の配列・pME008
  ├── logs/runtime.tsv           全ステップの所要時間（5種 × 11検体）
  └── refqc/                     小さい表（下記）

/mnt/scratch/tanoue/261016_Xiaoxi_eptAoligo_2wFeces_Metatranscriptome/   ← ROOT（中間ファイル）
  ├── ani/ pangenome/ mappability/ ref/ logs/ bakta/
  ├── trim/                      fastp 出力（11検体 × 2）
  ├── bam/                       bowtie2 の BAM 11本 + index
  ├── featureCounts/             検体ごとの .rds（11 × 2設定）
  ├── salmon_quant/              quant.sf ★ステージ5 で tximport が読む
  └── counts/                    featureCounts_q{0,10}.csv / feature_info.csv /
                                 salmon_{default,meta}_{reads,tpm,efflen}.csv / salmon_all.rds

/mnt/scratch/tanoue/genome/261016_Xiaoxi_eptAoligo_5strains/
  ├── fna/<株ID>.fna             444 で読み取り専用。md5 検証済み
  └── provenance.tsv             accession・md5・由来
```

**ステージ3 が project 側に置いた表（`refqc/`）**

| ファイル | 内容 |
|---|---|
| `eptA_cds_coordinates.tsv` | eptA 変異の CDS 内座標 |
| `stage3_bio.csv` | biotype 別の内訳 |
| `stage3_strain_cds.csv` | 株ごとの CDS カウントとシェア（除外基準の判定） |
| `stage3_vpi_st8.csv` | `Bth_VPI : Bth_St8` の比 |
| `stage3_assign_sum.csv` | featureCounts の割当内訳 |
| `salmon_mapping.tsv` | Salmon のマッピング内訳（CDS率 / decoy行き率 / ISR率） |
| `salmon_strain_share.tsv` | Salmon の株別シェア（2設定） |
| `salmon_vs_featurecounts.tsv` | 両者の一致度 |

**ゲノムとアノテーションの所在（`/db`。すべて取得済み。ダウンロード不要）**

| 株ID | accession | 由来 | ディレクトリ |
|---|---|---|---|
| `Bth_VPI` | GCF_000011065.1 | RefSeq | `/db/refseq/250328_ncbi_genome/Bacteroides_thetaiotaomicron_VPI-5482/GCF_000011065.1_ASM1106v1/` |
| `Bth_St8` | GCA_022835035.1 | DDBJ (BQNN01) | `/db/lab_genomes/05_CE91-68st/CE91_08_Bacteroides_thetaiotaomicron/GCA_022835035.1_ASM2283503v1/` |
| `Bun_8492` | GCF_000154205.1 | RefSeq | `/db/refseq/250328_ncbi_genome/Bacteroides_uniformis_ATCC_8492/GCF_000154205.1_ASM15420v1/` |
| `Odo_St23` | GCA_022835375.1 | DDBJ (BQOE01) | `/db/lab_genomes/05_CE91-68st/CE91_23_Odoribacteraceae_bacterium/GCA_022835375.1_ASM2283537v1/` |
| ~~`Odo_St24`~~ | GCA_022835415.1 | DDBJ (BQOG01) | **使わない**（§3） |

各ディレクトリに `_genomic.fna.gz` / `_genomic.gff.gz` / `_translated_cds.faa.gz` /
`_protein.faa.gz` / `_assembly_report.txt` が揃っている。

**既存 GFF の feature 数（2026-10-20 確認。Bakta 不要の根拠）**

| 株 | 由来 | CDS | rRNA | tRNA | tmRNA |
|---|---|---|---|---|---|
| Bth_VPI | RefSeq 2024-09 | 4,802 | 15 | 70 | 1 |
| Bth_St8 | DDBJ | 4,795 | 15 | 71 | 1 |
| Bun_8492 | RefSeq 2024-04 | 3,795 | 12 | 56 | 1 |
| Odo_St23 | DDBJ | 4,618 | 15 | 55 | 1 |

**4株すべてに rRNA・tRNA・tmRNA が入っている。Bakta の DB ダウンロードは不要。**

マウスゲノム（decoy）: `/db/gencode_mouse/GRCm39_M36`

---

## 付録 A. 前版から引き継いだ項目

前版 `CLAUDE_md_261016_previous.md` にあり、本版の本文から落ちている内容。
本プロジェクトでは既に実験が終わっているため適用できないものを含むが、**次の oligo 実験の
設計時に参照する**ため残す。

### A-1. 実験計画側（サンプリング前に決まってしまう事項）

1. **同一ペレットから DNA と RNA を分取**（AllPrep 等）。RNA/DNA 比が使えると正規化の悩みが大半消える
2. **spike-in を便重量あたり一定量で、抽出前に添加**（ERCC または外来菌体）。
   絶対量に変換でき、compositionality の問題が原理的に消える
3. **採便後ただちに凍結**。可能なら嫌気チャンバー内で採取。採取時刻を揃える（概日リズム）
4. **ケージを群あたり複数に分散し、ケージIDを記録**。gnotobiotic ではケージが定着の単位
   （※本プロジェクトは別アイソレータ飼育で群と交絡。§1 参照）
5. **深度**: rRNA残存50%を前提に、株レベルまで見るなら mRNA換算で 40〜60M read/sample
6. **株特異的 qPCR を使うなら、プライマーの交差反応を事前確認**。
   過去プロジェクトでは33株中11株のプライマーが複数株を同時に検出していた
7. **近縁株を同時に入れない。** 必要でなければ避け、必要なら事前に ANI と mappability を計算する。
   本プロジェクトの Odo_St23/St24 がまさにこの問題に当たった（§3）

### A-2. DESeq2 (A)/(B) の食い違いの実例

過去データでは (A) と (B) で結論が食い違う経路が複数あった。
**PTS は (A) で padj 1e-12、(B) では n.s.** — *Enterococcus* が3倍に増えただけだった。
(A) の結果を「発現制御の変化」と読んではいけない根拠。

### A-3. SortMeRNA をどうしても使う場合

SILVA ではなく、**本プロジェクトの株自身の rRNA 配列**（既存 GFF の rRNA feature から抽出）を
DB にし、`--fastx --aligned` で rRNA 側も残してカウントする。
ただし §6-6 のとおり、本プロジェクトでは使わない。
