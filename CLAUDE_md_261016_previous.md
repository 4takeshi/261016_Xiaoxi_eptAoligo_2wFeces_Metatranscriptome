# oligo5 metatranscriptome — 解析設計書

このファイルはプロジェクトルートに置く。Claude Code はこれを自動で読む。
**ここに書かれた判断は既に検討済み。再検討せず、そのまま従うこと。**
変更が必要だと判断した場合は、実行せずに理由を述べて確認を取る。

## 0. 実験と目的

- gnotobiotic mouse、5株 oligo-colonization、糞便 bacterial metatranscriptome
- 群間比較（食餌等）。各群 n≥5、複数ケージに分散、ケージIDを記録する
- 構成株:

| ID | 株 | 備考 |
|----|----|----|
| Bth_VPI | *B. thetaiotaomicron* VPI-5482 | type strain。プラスミド p5482 を忘れず含める |
| Bth_St8 | *B. thetaiotaomicron* CE-St8 | **Bth_VPI と同種。近縁** |
| Bun_8492 | *B. uniformis* ATCC 8492 | 別種。問題なし |
| Odo_St23 | *Odoribacter* CE-St23 | **Odo_St24 とほぼ同一ゲノム** |
| Odo_St24 | *Odoribacter* CE-St24 | **Odo_St23 とほぼ同一ゲノム** |

**5株中2ペアが近縁。これが本プロジェクト最大の技術的制約。**
ステージ1のゲートを通過するまで、ステージ2以降に進まないこと。

## 1. 絶対に守ること（過去プロジェクトで実際に起きた誤りの再発防止）

1. **`featureCounts` は必ず `-s 2`（reverse-stranded）。** Rsubread なら `strandSpecific=2`。
   根拠: 同系統の過去データで `ReadsPerGene.out.tab` の rev/fwd = 20〜55:1 を実測済み。
   未指定だとアンチセンスリードが約4%誤カウントされ、Ambiguity も増える。
2. **DESeq2 に入れるのは生カウントのみ。** カウントを係数で割ったり掛けたりして `round()` してから
   DESeq2 に渡してはならない。NB分散モデルの前提が壊れ、縮小した対象は検出力を失い、
   拡大した対象は偽陽性が増える。補正が必要なら `normalizationFactors()` に offset として渡す。
   過去データでは、この誤りで検定可能遺伝子の30%と DEG の40%を失っていた。
3. **CDS のみで DE 解析する。** GTF の `source == "Prodigal"` で絞る。
   rRNA/tRNA/ncRNA/CRISPR はカウントはするが QC 表に回す。
   過去データでは rRNA が assigned read の 37〜56% を占めたまま DESeq2 に入っていた。
4. **`enricher()` には必ず `universe` を渡す。** 背景は「検定された遺伝子（padj が NA でない）のうち
   対象集合に属するもの」。渡さないと「変動した機能」ではなく「その株が持っている機能」が出る。
5. **ゲノムの配列を連結して1本にしない。** ファイルを連結して1つの multi-FASTA にする。
   contig は独立したレコードのまま保つ。
6. **contig 名に株IDを入れる。** Bakta は contig を `contig_1` 等に改名するため、
   そのまま5株分を連結すると名前が衝突する。`03_build_reference.sh` が処理する。
7. **重複除去（Picard MarkDuplicates 等）をしない。** RNA-seq では高発現遺伝子が
   正当に重複リードを生むため、入れると高発現遺伝子だけが壊れる。
8. **SortMeRNA を使わない。** 全部マップしてから feature で CDS に絞る。
   rRNA% が QC 指標として自動的に得られる。どうしても使うなら SILVA ではなく
   Bakta が同定した本プロジェクトの株自身の rRNA 配列を DB にし、`--fastx --aligned` で
   rRNA 側も残してカウントする。

## 2. DESeq2 は2本立てる（1本では両方の問いに答えられない）

```
(A) 全株プール DESeq2    → 群集としての機能出力。菌量変化を含む。宿主が受け取る量
(B) 株ごと独立 DESeq2    → 株内の発現制御。菌量変化を含まない。各菌が何を変えたか
```

**(B) が主、(A) が従。** 過去データでは、両者で結論が食い違う経路が複数あった
（PTS は (A) で padj 1e-12、(B) では n.s. — *Enterococcus* が3倍に増えただけだった）。

(A) 単独でやると、低存在量株の dispersion が高存在量株から借りられ、偽陽性が出る。
(B) 単独でやると、群集レベルの機能出力が見えない。

**(B) は株ごとに size factor を持つため、株全体のベースラインずれが構造的に起きない。**
(A) に per-strain offset を足す方式は、株全体が丸ごとシフトして全経路が汚染されるので採らない。

付随事項:
- `lfcShrink(type="apeglm")` を入れる。FCでランク付けする箇所は必ず収縮後の値を使う
- prefiltering: 最小群サイズ以上のサンプルで count >= 10
- 検定は Wald（2群なので LRT と実質同等）
- design に cage を入れられるよう coldata に cage 列を持たせる

## 3. ステージとゲート

### ステージ0 — 環境確認
`00_check_env.sh`。足りないものを `mamba` で入れる。

### ステージ1 — リファレンスQC（**ここがゲート**）
`01_reference_qc.sh`, `02_mappability.sh`

出す数字:
- 全ペアの ANI と**アライン割合**（ANI 単独では判断できない）
- 各株の一意配列量（kb）
- pangenome の遺伝子有無行列（1株にしかない遺伝子の数）
- CDS ごとの mappability（実リード長で実測）と、混同相手の株

**判定基準**（近縁ペアごとに適用）:

| 一意配列 / 1株にしかない遺伝子 | 判定 |
|---|---|
| 100 kb 以上 / 70遺伝子以上 | 同一マウスで可。アクセサリ遺伝子で株差を論じる |
| 10〜100 kb | 境界。論じられる遺伝子を明示したうえで進む |
| 10 kb 未満 | **同一マウスに入れる意味がない。片方を外すか群を分ける** |

**この判定はヒトが行う。スクリプトは数字を出すところまで。**

### ステージ2 — リファレンス構築
`03_build_reference.sh`
CheckM2（完全性・コンタミ） → Bakta → eggNOG-mapper → contig改名 → 統合FASTA
→ **マウスゲノムを decoy として追加** → 統合GTF（株ID・CDSフラグ・KO付き）

マウス decoy は必須。過去データで `unmapped: too short` が 17〜59% あり、
宿主RNA が細菌遺伝子に誤マップしている可能性が高い。

### ステージ3 — サンプル処理
`04_reads.sh`
FastQC/MultiQC → fastp（`--trim_poly_g`, `--length_required 50`）
→ Bowtie2（または STAR `--alignIntronMax 1`）
→ samtools sort/index → featureCounts（`-s 2 -p --countReadPairs -B --fraction -M`）
→ **Salmon を並走**（selective alignment、decoy に全ゲノム）

近縁株があるため Salmon は必須。featureCounts の fraction mode は曖昧リードを 1/NH で
機械的に等分するだけだが、Salmon の EM は一意リードの情報を使って確率的に再配分する。
両者の相関を見て、乖離の大きい遺伝子にフラグを立てる。

### ステージ4 — QCと定着判定（**第2のゲート**）
`05_qc_table.R`

表に出す: サンプル別 mapping率 / rRNA% / 株別CDSリード数・シェア・unique率 / multimap率

**除外基準: 総CDSリード 200,000 未満、または全CDSリードの 0.1% 未満の株は
株レベル解析から外す。** 過去プロジェクトでは全8検体で1リードしかない株が
最後まで解析に残り、論文の主張を支える「4株」の一員になっていた。

### ステージ5 — 統計
`06_deseq2.R`。上記 §2 のとおり (A)(B) 両方。

### ステージ6 — 機能解析
`07_enrichment.R`
- **fgsea が主**（DESeq2 の Wald 統計量でランク付け）。ORA は従
- ORA には universe 必須
- アノテーションは KEGG pathway だけにしない。**KEGG module (ko00002)、CAZy/dbCAN（+ Bacteroides の PUL）、
  GMM/GBM、TCDB** を併用。*Bacteroides* 主体の群集なので多糖利用が主役になる可能性が高い
- **除外**: ko01100 / 01110 / 01120 / 01200 / 01210 / 01212 / 01230 等の global map。
  eggNOG が細菌に誤付与する真核・疾患 pathway（cytochrome P450、retinol、各種 signaling）
- **KEGG のラベルを鵜呑みにしない**。例: ko00195「Photosynthesis」の実体は F0F1 ATP synthase、
  ko00120/00121「bile acid biosynthesis」の実体は choloylglycine hydrolase（BSH）。
  図には実体を書く
- 有意経路はゲノム座標順に LFC をプロットし、1オペロンが効いているだけでないか確認する

## 4. 実験計画側（ドライではなく実験で決まる事項）

解析開始後では取り返せないので、サンプリング前に確定すること。

1. **同一ペレットから DNA と RNA を分取**（AllPrep 等）。RNA/DNA 比が使えると正規化の悩みが大半消える
2. **spike-in を便重量あたり一定量で、抽出前に添加**（ERCC または外来菌体）。
   絶対量に変換でき、compositionality の問題が原理的に消える
3. **採便後ただちに凍結**。可能なら嫌気チャンバー内で採取。採取時刻を揃える（概日リズム）
4. **ケージを群あたり複数に分散し、ケージIDを記録**。gnotobiotic ではケージが定着の単位
5. **深度**: rRNA残存50%を前提に、株レベルまで見るなら mRNA換算で 40〜60M read/sample
6. **株特異的 qPCR を使うなら、プライマーの交差反応を事前確認**。
   過去プロジェクトでは33株中11株のプライマーが複数株を同時に検出していた

## 5. ディレクトリ

```
oligo5_metatx/
├── CLAUDE.md
├── genomes/          入力: 5株の .fna（株IDをファイル名に）
├── bakta/            ステージ2出力
├── ref/              統合FASTA / index / GTF
├── ani/              ステージ1: ANI・一意配列
├── pangenome/        ステージ1: Panaroo
├── mappability/      ステージ1: CDSごとの区別可能性
├── fastq/            raw / trimmed
├── bam/
├── counts/
├── results/
└── logs/
```

## 6. 報告の作法

各ステージの終わりに、次を簡潔に報告すること。
- 出した数字（表で）
- 判断が必要な箇所と、その判断に必要な情報
- 想定と違った点

**ゲート（ステージ1終了時、ステージ4終了時）では必ず止まり、確認を取ってから次に進む。**
