# Genotyping Primer Designer for Deletion（欠失ジェノタイピング用プライマー設計）

*[English README](README.md)*

「どこからどこまでを欠失させるか」を指定すると、編集クローンと未編集クローンを
**バンドの長さの違い**で判定する PCR を設計します。プライマーは修復テンプレート
（ドナー）のホモロジーアームの外側に置かれ、WT 側・欠失側どちらの増幅産物長も
配列そのものから実測して報告します。

```
 ... [ 探索範囲 ][ arm ][      欠失させる領域      ][ arm ][ 探索範囲 ] ...
       ^ forward primer                                 reverse primer ^
                |<---- clearance（余裕） ---->|

        未編集アレル  ──────────────────────────────────────────  1,427 bp
        編集アレル    ───────────────/\───────────────────────     926 bp
                                     欠失（＋任意の挿入配列）
```

エクソンを抜く前提は一切ありません。欠失させるブロックは指定した任意の区間
（エクソン、プロモーター、調節領域、プラスミド上のスタッファー配列など）であり、
入力は Ensembl から取得したゲノム領域でも、手元の FASTA でも、SnapGene の
環状コンストラクトでも構いません。

## なぜプライマーの位置に制約があるのか

汎用のプライマー設計ツールでは扱えない部分であり、本ツールを分離した理由でもあります。

**どちらのプライマーもホモロジーアームに重なってはいけません。** ssODN ドナーや
ターゲティングベクターは切断点の直近配列をそのまま持っています。アームの内側に
アニールするプライマーは、ゲノムと同じようにドナーにもアニールします。その結果、
トランスフェクション試薬の残存が増幅され、未編集クローンまで「編集済み」と読めて
しまいます。`--arm` はそのアーム長で、欠失領域の両側それぞれこの塩基数を
プライマー探索から除外します。ドナーを使わない場合は `--arm 0` を指定します。

**アームの外側でも、切断点からの距離には意味があります。** Cas9 切断部位には
インデルが残り、想定以上に削れることもあります。アームの 15 bp 外側のプライマー
配列は接合部と一緒に壊れ得ますが、200 bp 外側ならまず壊れません。そのため
`--min-clearance`（既定 20 bp）より切断点に近いプライマーも**候補としては提示**
しますが（余裕のない領域ではそれしか選べないこともあるため）、必ず理由付きの
警告を添えます。黙って捨てることも、リスクを黙って通すこともしません。

順位付けも同じ優先順です。ユニークでないプライマー、指定した clearance の内側に
入るプライマーは、条件を満たすペアより必ず下に落とします。条件を満たすペアどうしは
Tm と GC、次に接合部からの距離、最後に増幅産物長で並べます。

## インストール

Python 3.10 以降。依存パッケージは numpy のみです。

```bash
git clone https://github.com/takeisue1270/Genotyping_Primer_Designer_for_Deletion.git
cd Genotyping_Primer_Designer_for_Deletion
pip install -e .
```

インストールせずにソースから直接実行する場合:

```bash
PYTHONPATH=src python -m genotyping_primers.cli --help
```

## 使い方

### 手元の配列ファイルを使う

```bash
genotyping-primers --sequence examples/example_locus.fa --delete 5501-6000
```

```
example_locus   examples/example_locus.fa
  12,000 bp linear   coordinates: example_locus (1-based)

  deletion   5,501-6,000   500 bp removed
  net        the edited allele is 500 bp shorter
  excluded   43 bp homology arm on each side; 20 bp of clearance asked for beyond it

  Recommended pair
    example_locus_F    ACGACGTAAGGTGAACTTTCGGTCT
                       25 nt  Tm 60.0 C  GC 48 %  + strand  5,175-5,199  258 bp clear of the edit
                       local: 1 exact, 1 within 2 mismatch(es)
    example_locus_R    GCCACCCAGATTCATCCACTTAATGC
                       26 nt  Tm 60.0 C  GC 50 %  - strand  6,595-6,620  551 bp clear of the edit
                       local: 1 exact, 1 within 2 mismatch(es)
    bands: unedited 1,446 bp / edited 946 bp   (difference 500 bp)
```

### 参照ゲノムの領域を使う

生物種と、ゲノム座標での欠失範囲を指定します。周辺配列は Ensembl から自動取得され、
プライマーの特異性は [GGGenome](https://gggenome.dbcls.jp/) でそのアセンブリ全体に対して
チェックされます。

```bash
genotyping-primers --species mouse --delete 11:69,000,000-69,001,000 -o designs
```

`mouse`、`human`、`rat`、`zebrafish`、`fly`、`worm`、`yeast` などの通称と、
Ensembl の正式名（`mus_musculus` など）のどちらでも指定できます。

### プラスミドの場合

コンストラクトは環状で、かつ問われているのは「ゲノム全体でユニークか」ではなく
「**そのプラスミドの中でユニークか**」です。`--specificity local`（生物種を
指定しない場合の既定値）はまさにその問いに答え、ネットワークには一切アクセス
しません。複製起点をまたぐ欠失も指定できます。

```bash
genotyping-primers --sequence pMyVector.gb --circular --delete 4950-250
```

```
  deletion   4,950-250   501 bp removed
  Recommended pair
    ...  + strand  4,574-4,597  309 bp clear of the edit
    ...  - strand    775-800    481 bp clear of the edit
    bands: unedited 1,427 bp / edited 926 bp   (difference 501 bp)
```

環状かどうかは GenBank の `LOCUS` 行から読み取るため、`--circular` が必要なのは
FASTA を渡すときだけです。

### 欠失部に配列を残す場合

`--insert` は欠失領域の代わりに残す配列（タグ、制限酵素サイト、バーコードなど）を
指定し、バンド長の計算にも反映されます。2 本のバンドの差は「欠失長 − 挿入長」
となり、これがレポートの net の値です。

```bash
genotyping-primers --sequence locus.fa --delete 5501-6000 --insert GGATCCAAGCTTGA
#   net   the edited allele is 486 bp shorter   （500 bp 欠失・14 bp 挿入）
```

## 出力

`-o` を付けなければレポートを表示するだけで、ファイルは作られません。付けた場合:

```
designs/
├── example_locus_genotyping_primers.tsv   プライマー 1 本につき 1 行（発注用）
├── example_locus_unedited.gb              未編集アレル（プライマーとアームを注釈）
└── example_locus_edited.gb                編集後アレル
```

TSV には配列、長さ、Tm、GC、位置、clearance、特異性チェック結果、両方のバンド長、
すべての警告が入ります。GenBank 2 ファイルは SnapGene でそのまま開けます。
未編集側のマップには欠失領域・両側のホモロジーアーム・プライマー結合部位が、
編集側のマップには検出対象アレル上の同じプライマーが描かれます。

レポートと TSV の座標はすべて**入力時の座標系**です。ゲノム領域を取得した場合は
ゲノム座標、ファイルを読んだ場合はファイル内座標、プラスミドの場合は元の番号のまま
報告します（環状配列は探索の余裕を確保するため内部で回転させていますが、表示時に
戻します）。例外は GenBank マップで、回転させたプラスミドは回転後の座標で書き出し、
`DEFINITION` 行に「元の何番目の塩基が 1 番になったか」を明記します。

## 警告の読み方

各警告はそれが該当するペアに付きます。ペアではなく「編集そのもの」の性質による
ものは **Notes** にまとめて出ます。

| 警告 | 意味 | 対応 |
| --- | --- | --- |
| *anneals N bp from the homology arm* | `--min-clearance` より接合部に近い。ドナーには重ならないが、接合部のインデルが届き得る。 | `--span` を広げて候補を増やす。納得の上で使うのも可。 |
| *no homology arm was declared* | `--arm 0` のため欠失領域そのものしか避けていない。 | 修復テンプレートを使うなら `--arm` にアーム長を指定して再実行。 |
| *specificity is unknown* | チェックが実行されなかった、またはサービスが応答しなかった。 | 再実行するか手動で確認。未チェックのプライマーを「問題なし」とは決して報告しません。 |
| *is not unique* | 許容ミスマッチ内に他のゲノム部位がある。 | 別のペアを使うか `--span` を広げる。 |
| *3' N nt occur N times* | 全長ではユニークだが 3' 末端が反復配列。ミスプライミングはここから起こる。 | 代替ペアを優先。 |
| *the two bands differ by only N bp* | 差が `--min-band-difference`（既定 50 bp）未満。通常のアガロースでは分離しない。 | 欠失範囲を広げる、高濃度ゲルを使う、接合部をシークエンスする。 |
| *the edit does not change the length* | 挿入配列と欠失長が同じで、サイズ差による判定が原理的に不可能。 | 接合部のシークエンス、または挿入配列内部にプライマーを置く。 |
| *the unedited band is N bp* | 3 kb 超。反応失敗と欠失アレルの区別がつきにくい。 | ポジティブコントロールを併走させるか `--span` を下げる。 |
| *the edited band is only N bp* | 100 bp 未満。ゲルから流れ出やすい。 | `--span` を上げてプライマーを外側へ。 |

## 特異性チェック

`--specificity` で「このプライマーはユニークか」の確認方法を選びます。

| バックエンド | 対象範囲 | ミスマッチ部位も検出 | ネットワーク |
| --- | --- | --- | --- |
| `gggenome` | アセンブリ全体 | する（`--mismatches`、既定 2） | 必要 |
| `fasta` | ローカルのゲノム FASTA | しない（完全一致のみ＋3' 末端 15 nt の出現数を別途カウント） | 不要 |
| `local` | 与えた配列の中だけ | する | 不要 |
| `none` | チェックしない | — | 不要 |
| `auto`（既定） | `--genome-fasta` があれば `fasta`、種に DB があれば `gggenome`、それ以外は `local` | | |

重要な点が 2 つあります。

* **ゲノムを対象とする設計で `local` を使った場合、「ユニークでない」ことは示せても
  「ユニークである」ことは示せません。** レポートにもその旨を明記し、チェック済みの
  ように見せることはしません。プラスミドでは逆に、これが妥当なチェックです。
* **GGGenome は存在しない DB 名を渡すとエラーではなくヒトゲノムの結果を返します。**
  本ツールは応答が期待した生物種のものかを毎回検証するため、`--gggenome-db` の
  打ち間違いは黙って hg38 で採点されるのではなく、明示的に失敗します。既定以外の
  アセンブリ（GRCm38 系統なら `mm10` など）を使いたい場合に `--gggenome-db` を指定します。

## オプション一覧

| オプション | 既定値 | |
| --- | --- | --- |
| `--sequence PATH` | | FASTA / GenBank(SnapGene) / 生配列 |
| `--record NAME` | 先頭 | マルチレコード FASTA のどのレコードを使うか |
| `--species NAME` | | Ensembl から取得し、そのアセンブリで特異性を確認 |
| `--region CHR:START-END` | | 取得する範囲を明示指定 |
| `--pad BP` | arm + span + 200 | `--delete` の周囲に取得する余白 |
| `--circular` | GenBank の記載に従う | 環状配列として扱う |
| `--delete SPEC` | **必須** | `4800-5300`、`4800..5300`、`4800+500`、`11:69,000,000-69,001,000` |
| `--coordinates` | auto | `--delete` の解釈：`local`（配列内座標）／`genomic`（ゲノム座標）／自動判定 |
| `--insert SEQ` | なし | 欠失部に残す配列 |
| `--arm BP` | 43 | ホモロジーアーム長。両側この塩基数を探索から除外 |
| `--min-clearance BP` | 20 | 除外域からこの距離未満なら警告 |
| `--span BP` | 1000 | 除外域の外側どこまで探索するか |
| `--tm C` | 60 | 目標 Tm（最近接塩基対法、SantaLucia 1998） |
| `--tm-range LO,HI` | 58,62 | 許容 Tm 範囲 |
| `--max-tm-diff C` | 3 | ペア内の Tm 差の上限 |
| `--length LO,HI` | 20,30 | プライマー長 |
| `--gc LO,HI` | 0.4,0.6 | GC 含量 |
| `--specificity` | auto | `auto` / `gggenome` / `fasta` / `local` / `none` |
| `--gggenome-db DB` | 種から自動 | GGGenome のデータベース名 |
| `--genome-fasta PATH` | | オフラインでゲノム全体を数えるためのローカル FASTA |
| `--mismatches N` | 2 | オフターゲット計数時に許容するミスマッチ数 |
| `-o, --outdir DIR` | 表示のみ | TSV とマップの出力先 |
| `--alternatives N` | 2 | 併記する次点ペアの数 |
| `--no-maps` | | GenBank を出力しない |
| `--min-band-difference BP` | 50 | これ未満を「分離不可」と判定 |
| `-v, --verbose` | | 進捗ログを表示 |

終了コード: `0` 設計成功、`1` 実行条件が整わない（座標不正、ファイルが読めない、
その種の DB がない）、`2` 実行は正常だがペアが見つからなかった。

上記のほか、同一塩基 5 連続、相手プライマーとの 3' 末端相補 5 塩基以上、
GC クランプ不適（末端 5 塩基中の G/C が 1〜3 個でない）、曖昧塩基との重なりを
持つ候補も除外されます。

## ライブラリとして使う

```python
from genotyping_primers import DesignConfig, design, prepare, resolve_span, target_from_file

target, notes = target_from_file("locus.fa")
target, deletion = prepare(target, *resolve_span(target, "5501-6000"))
result = design(target, deletion, DesignConfig(homology_arm=43, specificity="local"))

best = result.best
print(best.forward.sequence, best.reverse.sequence)
print(best.wt_product, best.ko_product, best.warnings)
```

`design` は「設計結果が良くなかった」だけでは例外を投げません。見つかったペアを
警告付きで返し、できなかったことは notes に記録します。例外を投げるのは、そもそも
実行条件が整わないときだけです。

## このツールがしないこと

プライマーを選ぶだけで、ガイド RNA の設計はしませんし、切断位置の良し悪しも
判断しません。予測しているのはアニーリングであって増幅ではないため、ここでの
条件をすべて満たすペアでも、GC リッチな鋳型や二次構造で失敗することはあります。
特異性チェックも「部位の数」を数えているだけで、ポリメラーゼが実際に伸長するか
どうかまでは分かりません。下流のコストが大きい実験では、最終的なペアを
in-silico PCR でも確認してください。

## テスト

```bash
pip install -e ".[dev]"
pytest
```

テストは完全にオフラインで動きます。固定シードから合成配列を作るため期待座標が
厳密に決まり、報告された増幅産物長は「その設計が増幅すると主張している配列」から
毎回測り直して検証しています（ツール自身の計算を信用しません）。座標変換
（ファイル内・ゲノム・環状回転）、GenBank の往復、FASTA スキャナのチャンク境界処理、
GGGenome クライアント（hg38 フォールバックを含むスタブ応答）もカバーしています。

## ライセンス

MIT。[LICENSE](LICENSE) を参照してください。
