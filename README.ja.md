# Genotyping Primer Designer for Deletion（欠失ジェノタイピング用プライマー設計）

*[English README](README.md)*

欠失させる区間を指定すると、その欠失を持つアレルと持たないアレルを増幅産物長で
判別する PCR を設計します。プライマーは修復テンプレートのホモロジーアームの外側に
配置します。増幅産物長は両アレルの配列から実測して報告します。

```
 ... [ 探索範囲 ][ arm ][      欠失させる領域      ][ arm ][ 探索範囲 ] ...
       ^ forward primer                                 reverse primer ^
                |<---- clearance（余裕） ---->|

        未編集アレル  ──────────────────────────────────────────  1,427 bp
        編集アレル    ───────────────/\───────────────────────     926 bp
                                     欠失（＋任意の挿入配列）
```

欠失させる区間は任意です。エクソン、プロモーター、調節領域、プラスミド上の
スタッファー配列のいずれでも指定できます。入力は Ensembl から取得したゲノム領域、
手元の FASTA、SnapGene の環状コンストラクトのいずれでも構いません。

## プライマーの位置に制約を設ける理由

汎用のプライマー設計ツールが扱わない点であり、本ツールを作成した理由です。

**第一に、どちらのプライマーもホモロジーアームに重ねてはいけません。** ssODN ドナーや
ターゲティングベクターは、切断点の直近配列をそのまま持っています。アームの内側に
アニールするプライマーは、ゲノムだけでなくドナーにもアニールします。その結果、
残存したトランスフェクション試薬が増幅され、未編集クローンも編集済みと判定されます。
`--arm` にアーム長を指定すると、欠失領域の両側のその範囲を探索から除外します。
ドナーを使わない場合は `--arm 0` とします。

**第二に、アームの外側でも切断点からの距離が問題になります。** Cas9 の切断部位には
インデルが残り、想定以上に削れることもあります。アームの 15 bp 外側のプライマー配列は
接合部とともに失われ得ますが、200 bp 外側であればまず失われません。そこで
`--min-clearance`（既定 20 bp）より切断点に近いプライマーも候補として提示しますが、
理由を付した警告を添えます。余裕のない領域ではそれしか選べないためです。

順位付けにも同じ優先順を用います。ユニークでないプライマーと、指定した clearance の
内側に入るプライマーは、条件を満たすペアより下位に置きます。条件を満たすペアどうしは、
Tm と GC、接合部からの距離、増幅産物長の順に評価します。

## インストール

Python 3.10 以降が必要です。依存パッケージは numpy のみです。

```bash
git clone https://github.com/takeisue1270/Genotyping_Primer_Designer_for_Deletion.git
cd Genotyping_Primer_Designer_for_Deletion
pip install -e .
```

インストールせずにソースから実行することもできます。

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

生物種と、ゲノム座標での欠失範囲を指定します。周辺配列は Ensembl から取得します。
プライマーの特異性は [GGGenome](https://gggenome.dbcls.jp/) を用いて、そのアセンブリ
全体に対して確認します。

```bash
genotyping-primers --species mouse --delete 11:69,000,000-69,001,000 -o designs
```

`mouse`、`human`、`rat`、`zebrafish`、`fly`、`worm`、`yeast` などの通称と、
Ensembl の正式名（`mus_musculus` など）のいずれでも指定できます。

### プラスミドの場合

プラスミドで問うべきことは、ゲノム全体での一意性ではなく、そのプラスミド内での
一意性です。`--specificity local`（生物種を指定しない場合の既定値）がこれに答えます。
ネットワークにはアクセスしません。複製起点をまたぐ欠失も指定できます。

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

環状かどうかは GenBank の `LOCUS` 行から読み取ります。`--circular` が必要なのは
FASTA を渡す場合のみです。

### 欠失部に配列を残す場合

`--insert` で、欠失領域の代わりに残す配列を指定します。タグ、制限酵素サイト、
バーコードなどが該当します。2 本のバンドの差は「欠失長 − 挿入長」であり、これが
レポートの net の値です。

```bash
genotyping-primers --sequence locus.fa --delete 5501-6000 --insert GGATCCAAGCTTGA
#   net   the edited allele is 486 bp shorter   （500 bp 欠失・14 bp 挿入）
```

## 出力

`-o` を付けない場合はレポートを表示するだけで、ファイルは作成しません。
付けた場合は次の 3 ファイルを出力します。

```
designs/
├── example_locus_genotyping_primers.tsv   プライマー 1 本につき 1 行（発注用）
├── example_locus_unedited.gb              未編集アレル（プライマーとアームを注釈）
└── example_locus_edited.gb                編集後アレル
```

TSV には、配列、長さ、Tm、GC、位置、clearance、特異性チェック結果、両バンド長、
すべての警告が入ります。GenBank 2 ファイルは SnapGene でそのまま開けます。
未編集側には欠失領域、両側のホモロジーアーム、プライマー結合部位を、編集側には
検出対象アレル上の同じプライマーを描画します。

レポートと TSV の座標は、いずれも入力時の座標系です。ゲノム領域を取得した場合は
ゲノム座標、ファイルを読んだ場合はファイル内座標、プラスミドの場合は元の番号で
報告します。環状配列は探索範囲を確保するため内部で回転させますが、表示時に戻します。
例外は GenBank マップです。回転させたプラスミドは回転後の座標で書き出し、元の何番目の
塩基が 1 番になったかを `DEFINITION` 行に明記します。

## 警告の読み方

警告は該当するペアに付きます。ペアではなく編集そのものの性質による警告は、
**Notes** にまとめます。

| 警告 | 意味 | 対応 |
| --- | --- | --- |
| *anneals N bp from the homology arm* | `--min-clearance` より接合部に近い。ドナーには重ならないが、接合部のインデルが届き得る。 | `--span` を広げて候補を増やす。了解の上で使うこともできる。 |
| *no homology arm was declared* | `--arm 0` のため、欠失領域のみを避けている。 | 修復テンプレートを使う場合は `--arm` にアーム長を指定して再実行する。 |
| *specificity is unknown* | チェックが実行されなかった、またはサービスが応答しなかった。 | 再実行するか手動で確認する。未チェックのプライマーを問題なしとは報告しない。 |
| *is not unique* | 許容ミスマッチ内に他のゲノム部位がある。 | 別のペアを使うか `--span` を広げる。 |
| *3' N nt occur N times* | 全長では一意だが、3' 末端が反復配列である。ミスプライミングはここから起こる。 | 代替ペアを優先する。 |
| *the two bands differ by only N bp* | 差が `--min-band-difference`（既定 50 bp）未満で、通常のアガロースでは分離しない。 | 欠失範囲を広げる、高濃度ゲルを用いる、接合部をシークエンスする。 |
| *the edit does not change the length* | 挿入長と欠失長が等しく、サイズ差による判定が原理的にできない。 | 接合部をシークエンスするか、挿入配列内部にプライマーを置く。 |
| *the unedited band is N bp* | 3 kb 超で、反応失敗と欠失アレルの区別がつきにくい。 | ポジティブコントロールを併走させるか `--span` を下げる。 |
| *the edited band is only N bp* | 100 bp 未満で、ゲルから流出しやすい。 | `--span` を上げてプライマーを外側へ移す。 |

## 特異性チェック

`--specificity` で、プライマーの一意性を確認する方法を選びます。

| バックエンド | 対象範囲 | ミスマッチ部位も検出 | ネットワーク |
| --- | --- | --- | --- |
| `gggenome` | アセンブリ全体 | する（`--mismatches`、既定 2） | 必要 |
| `fasta` | ローカルのゲノム FASTA | しない（完全一致のみ。3' 末端 15 nt の出現数は別途カウント） | 不要 |
| `local` | 与えた配列の中のみ | する | 不要 |
| `none` | チェックしない | — | 不要 |
| `auto`（既定） | `--genome-fasta` があれば `fasta`、種に DB があれば `gggenome`、それ以外は `local` | | |

注意すべき点が 2 つあります。

* ゲノムを対象とする設計で `local` を用いた場合、一意でないことは示せますが、
  一意であることは示せません。レポートにもその旨を明記します。プラスミドでは、
  これが妥当なチェックです。
* GGGenome は、存在しない DB 名に対してエラーではなくヒトゲノムの結果を返します。
  本ツールは応答が期待した生物種のものかを毎回検証するため、`--gggenome-db` の
  誤記は hg38 で採点されるのではなく、明示的に失敗します。既定以外のアセンブリ
  （GRCm38 系統の `mm10` など）を使う場合に `--gggenome-db` を指定します。

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
| `--arm BP` | 43 | ホモロジーアーム長。両側のこの範囲を探索から除外 |
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
| `--min-band-difference BP` | 50 | これ未満を分離不可と判定 |
| `-v, --verbose` | | 進捗ログを表示 |

終了コードは 3 種類です。`0` は設計成功、`1` は実行条件の不備（座標不正、ファイルが
読めない、その種の DB がない）、`2` は実行は正常だがペアが見つからなかった場合です。

上記のほか、同一塩基 5 連続、相手プライマーとの 3' 末端相補 5 塩基以上、GC クランプ
不適（末端 5 塩基中の G/C が 1〜3 個でない）、曖昧塩基との重なりを持つ候補も除外します。

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

`design` は、設計結果が良くないというだけでは例外を送出しません。見つかったペアを
警告付きで返し、できなかったことは notes に記録します。例外を送出するのは、実行条件が
整わない場合のみです。

## 本ツールが行わないこと

本ツールはプライマーを選ぶだけです。ガイド RNA の設計は行わず、切断位置の適否も
判断しません。予測しているのはアニーリングであって増幅ではないため、すべての条件を
満たすペアでも、GC リッチな鋳型や二次構造によって失敗することがあります。特異性
チェックも部位の数を数えているだけで、ポリメラーゼが実際に伸長するかどうかは
分かりません。下流のコストが大きい実験では、最終的なペアを in-silico PCR でも
確認してください。

## テスト

```bash
pip install -e ".[dev]"
pytest
```

テストはすべてオフラインで動作します。固定シードから合成配列を作るため期待座標が
厳密に定まります。報告された増幅産物長は、その設計が増幅すると主張している配列から
毎回測り直して検証しています。座標変換（ファイル内・ゲノム・環状回転）、GenBank の
往復、FASTA スキャナのチャンク境界処理、GGGenome クライアント（hg38 フォールバックを
含むスタブ応答）も対象としています。

## ライセンス

MIT です。[LICENSE](LICENSE) を参照してください。
