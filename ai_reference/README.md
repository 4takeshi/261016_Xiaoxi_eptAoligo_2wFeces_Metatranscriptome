# AI からの参考ファイル

解析の設計を検討する段階で、claude.ai が作成した notebook。
**実行していない参考資料**であり、本解析のパイプライン（`1_` 〜 `7_`）には含まれない。

| ファイル | 中身 | 本解析での扱い |
|---|---|---|
| `00_01_env_and_reference_qc.ipynb` | ステージ0–1（環境確認・ANI・一意配列・pangenome・mappability）の草案。パスは `/mnt/hdd1/tanoue/oligo5_metatx`、5株（`Odo_St24` を含む）想定 | `1_genome_copy.ipynb` と `2_reference_qc.ipynb` として書き直した。**4株構成に変更**し、`Odo_St24` は外した |
| `3_build_reference_claudeai_original.ipynb` | ステージ2（リファレンス構築）の草案 | `3_build_reference.ipynb` として書き直した |

★ 中のパス・株構成・パラメータは**現行の解析と一致しない**。設計の経緯を残すためだけに置いてある。
