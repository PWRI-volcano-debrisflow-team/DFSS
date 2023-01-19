---
author: 清水武志
affiliation: 土木研究所　火山・土石流チーム
date: 2023/1/19
---


# はじめに

古恵川のデータ（`../03_furuekawa` のこと。以下、オリジナルコードという。）を参考に、計算に必要な入出力ファイルをわかりやすくするためにフォルダ構成などを再整理したデータである。

以下にオリジナルコードとの相違点を示す。
ただし、前処理である地形データのGIS処理および後処理である可視化に関するプログラムは対象外とした。


## ソースコードの変更点について

`./bin` に fortran のソースコードを収める。オリジナルコードとそのパッチファイルを収めた。

下流域地形の前処理プログラムである `01_mk_watershedConfigurationFiles.f90` および `02_mk_floodplainConfigurationFiles.f90` は、
上流域地形の前処理プログラムと区別しやすくするために、ファイル名冒頭の数値を10番台に変更した。

RR,DR,DFについては、起動直後に `openMP` で使用するコア数確認コードの追加、入出力ファイルのパス名を変更した。
変更点の詳細は `patch` ファイルを参照してほしい。
オリジナルコードから、以下のようにすれば、変更後のファイルが得られる。
```
patch -u < RR.path
```
これにより、`RR_ver1.0.f90` の内容が `RR_ver1.01.f90` と同様に変更される。
元に戻すときは、
```
patch -R < RR.path
```
とする。

ソースコードは 例えば、`gfortran` を使用する場合、以下のようにコンパイルして使用できる。
```
gfortran 01_mk_topographyFiles.f90 -o 01_mk_topographyFiles.exe
gfortran 02_mk_watershedConfigurationFiles.f90 -o 02_mk_watershedConfigurationFiles.exe
gfortran 03_mk_streamConfigurationFiles.f90 -o 03_mk_streamConfigurationFiles.exe
gfortran 11_mk_floodplainConfigurationFiles.f90 -o 11_mk_floodplainConfigurationFiles.exe
gfortran 12_mk_streamFloodplainConnectionFiles.f90 -o 12_mk_streamFloodplainConnectionFiles.exe
gfortran -fopenmp RR_ver_1.01.f90 -o RR_ver_1.01_omp.exe
gfortran -fopenmp DR_ver_1.01.f90 -o DR_ver_1.01_omp.exe
gfortran -fopenmp DF_ver_1.01.f90 -o DF_ver_1.01_omp.exe
```


## サンプルデータについて

`./*.asc, ./*.txt` が DFSS の入力ファイルである。

./input フォルダはこのサンプルを動かすためには使用しないが、上の `./*.asc` 等を入力ファイルを作成するために、
ユーザが外部から収集・準備する以下のデータの保存を想定したフォルダである。

+ 流出解析を行う流域の下流末端の座標 ( `./input/coordinate.txt or ./input/outletPoint_furue.kml` ）
+ 国土地理院数値標高データ（10m B） `FG-GML-XXXX-XX-DEM10B.zip or .xml`
+ XRAINによる雨量データ
+ 下流域の計算範囲指定ファイル ( `./input/area_ll.kml` ）

`./input` フォルダにおける地形データと雨量データは利用規約に従って保存していない。

以下のように処理を行なうことで、DFSS の計算に必要な入力ファイルが偉られる。

+ 地形データは、[Open Source GRASS GIS](https://grass.osgeo.org/) によるGIS処理による測地系変換やリサンプリング、傾斜処理を実施して、`./dem_10m.asc, ./dir_10m.asc, ./dem_30m.asc, ./dir_30m.asc` に変換
+ 雨量データは、DFSS形式の雨量データに `../rain_Aso.txt` に変換

作成方法の詳細は `../03_furuegawa/01_SET_BASIN` および `../03_furuegawa/02_MK_TOPO` を参照すること。なお、傾斜データ (`./dir_*m.asc`) で `GRASS GIS` による `r.watershed, r.fill.dir` によるアルゴリズムで得ることを想定している。他の GIS ソフトウェアの処理では上手く動作しない例が確認されている。

さらに、計算にはパラメータファイル（`./topographyConfiguration.txt, ./RR_input.txt` 等）が必要である。



## 計算結果保存フォルダについて

`bin` のコードはフォルダの有無を判定して自動生成は行わない。

そのため、計算前に `./output_RR, ./output_RR_misc, ./output_DR, ./output_DF` をあらかじめ作成しなければならない。



## 並列計算で使用するコア数の設定と実行方法について

計算を実行する前に `openMP` による並列計算で使用する`CPU`のコア数をwindows のコマンドプロンプト（`cmd.exe`）や macOS や Linux のターミナル (`zsh`） によって、環境変数と設定する必要がある。

`./` がカレントディレクトリの場合に、コアを8つ使用するには以下のように実行すればよい。
```
set OMP_NUM_THREADS=8
bin/01_mk_topographyFiles.exe
bin/02_mk_watershedConfigurationFiles.exe
bin/03_mk_streamConfigurationFiles.exe
bin/11_mk_floodplainConfigurationFiles.exe
bin/12_mk_streamFloodplainConnectionFiles.exe
bin\RR_ver_1.1_omp.exe
bin/DR_ver_1.1_omp.exe
bin/DF_ver_1.1_omp.exe
```

なお、`cmd.exe` の場合、パスの区切り文字を `./` が不要で、`/` は `\` に変更する必要があるかもしれない（バッチファイルにまとめる場合は `\\` ）。